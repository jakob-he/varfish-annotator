package com.github.bihealth.varfish_annotator.init_db;

import com.github.bihealth.varfish_annotator.VarfishAnnotatorException;
import com.github.bihealth.varfish_annotator.utils.VariantDescription;
import com.github.bihealth.varfish_annotator.utils.VariantNormalizer;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import java.io.File;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.SQLException;

/** Base class for gnomAD import. */
public final class CADDImporter {

  /** The name of the table in the database. */
  public static final String TABLE_NAME = "CADD";

  /** The JDBC connection. */
  private final Connection conn;

  /** Path to CADD tabix path. */
  private final String caddTabixpath;

  /** Helper to use for variant normalization. */
  private final String refFastaPath;

  /**
   * Construct the <tt>CADDImporter</tt> object.
   *
   * @param conn Connection to database
   * @param caddTabixpath Path to CADD tabix file.
   */
  public CADDImporter(Connection conn, String caddTabixpath, String refFastaPath) {
    this.conn = conn;
    this.caddTabixpath = caddTabixpath;
    this.refFastaPath = refFastaPath;
  }

  /** Execute CADD import. */
  public void run() throws VarfishAnnotatorException {
    System.err.println("Re-creating table in database...");
    recreateTable();

    System.err.println("Importing CADD...");
    final VariantNormalizer normalizer = new VariantNormalizer(refFastaPath);
    String prevChr = null;
    try (VCFFileReader reader = new VCFFileReader(new File(caddTabixpath), true)) {
      final CloseableIterator<VariantContext> it;
      it = reader.iterator();

      while (it.hasNext()) {
        final VariantContext ctx = it.next();
        if (!ctx.getContig().equals(prevChr)) {
          System.err.println("Now on chrom " + ctx.getContig());
        }
        importVariantContext(normalizer, ctx);
        prevChr = ctx.getContig();
      }
    } catch (SQLException e) {
      throw new VarfishAnnotatorException(
          "Problem with inserting into " + TABLE_NAME + " table", e);
    }

    System.err.println("Done with importing CADD...");
  }

  /**
   * Re-create the gnomAD table in the database.
   *
   * <p>After calling this method, the table has been created and is empty.
   */
  private void recreateTable() throws VarfishAnnotatorException {
    final String dropQuery = "DROP TABLE IF EXISTS " + TABLE_NAME;
    try (PreparedStatement stmt = conn.prepareStatement(dropQuery)) {
      stmt.executeUpdate();
    } catch (SQLException e) {
      throw new VarfishAnnotatorException("Problem with DROP TABLE statement", e);
    }

    final String createQuery =
        "CREATE TABLE "
            + TABLE_NAME
            + "("
            + "release VARCHAR(10) NOT NULL, "
            + "chrom VARCHAR(20) NOT NULL, "
            + "start INTEGER NOT NULL, "
            + "end INTEGER NOT NULL, "
            + "ref VARCHAR("
            + InitDb.VARCHAR_LEN
            + ") NOT NULL, "
            + "alt VARCHAR("
            + InitDb.VARCHAR_LEN
            + ") NOT NULL, "
            + "RawScore DOUBLE NOT NULL,"
            + "PHRED DOUBLE NOT NULL"
            + ")";
    try (PreparedStatement stmt = conn.prepareStatement(createQuery)) {
      stmt.executeUpdate();
    } catch (SQLException e) {
      throw new VarfishAnnotatorException("Problem with CREATE TABLE statement", e);
    }

    final ImmutableList<String> indexQueries =
        ImmutableList.of(
            "CREATE PRIMARY KEY ON " + TABLE_NAME + " (release, chrom, start, ref, alt)",
            "CREATE INDEX ON " + TABLE_NAME + " (release, chrom, start, end)");
    for (String query : indexQueries) {
      try (PreparedStatement stmt = conn.prepareStatement(query)) {
        stmt.executeUpdate();
      } catch (SQLException e) {
        throw new VarfishAnnotatorException("Problem with CREATE INDEX statement", e);
      }
    }
  }

  /** Insert the data from <tt>ctx</tt> into the database. */
  @SuppressWarnings("unchecked")
  private void importVariantContext(VariantNormalizer normalizer, VariantContext ctx)
      throws SQLException {
    final String insertQuery =
        "MERGE INTO "
            + TABLE_NAME
            + " (release, chrom, start, end, ref, alt, RawScore, PHRED)"
            + "VALUES ('GRCh37', ?, ?, ?, ?, ?, ?, ?)";

    final int numAlleles = ctx.getAlleles().size();
    for (int i = 1; i < numAlleles; ++i) {
      final VariantDescription rawVariant =
          new VariantDescription(
              ctx.getContig(),
              ctx.getStart() - 1,
              ctx.getReference().getBaseString(),
              ctx.getAlleles().get(i).getBaseString());
      final VariantDescription finalVariant = normalizer.normalizeInsertion(rawVariant);

      final PreparedStatement stmt = conn.prepareStatement(insertQuery);
      stmt.setString(1, finalVariant.getChrom());
      stmt.setInt(2, finalVariant.getPos() + 1);
      stmt.setInt(3, finalVariant.getPos() + finalVariant.getRef().length());
      stmt.setString(4, finalVariant.getRef());
      stmt.setString(5, finalVariant.getAlt());
      stmt.setDouble(6, ctx.getCommonInfo().getAttributeAsDouble("RawScore", 0.0));
      stmt.setDouble(7, ctx.getCommonInfo().getAttributeAsDouble("PHRED", 0.0));

      stmt.executeUpdate();
      stmt.close();
    }
  }
}
