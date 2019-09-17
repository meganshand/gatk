package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.BufferedLineReader;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.Metadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.MetadataUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.config.ConfigFactory;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Simple data structure to pass and read/write a List of {@link SimpleCount} objects.
 * Supports both TSV and HDF5.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SimpleCountCollection extends AbstractSampleLocatableCollection<SimpleCount> {
    private static final int DEFAULT_FEATURE_QUERY_LOOKAHEAD_IN_BP = 1_000_000;

    //note to developers: repeat the column headers in Javadoc so that they are viewable when linked
    /**
     * CONTIG, START, END, COUNT
     */
    enum SimpleCountTableColumn {
        CONTIG,
        START,
        END,
        COUNT;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }
    
    private static final Function<DataLine, SimpleCount> SIMPLE_COUNT_RECORD_FROM_DATA_LINE_DECODER = dataLine -> {
        final String contig = dataLine.get(SimpleCountTableColumn.CONTIG);
        final int start = dataLine.getInt(SimpleCountTableColumn.START);
        final int end = dataLine.getInt(SimpleCountTableColumn.END);
        final int count = dataLine.getInt(SimpleCountTableColumn.COUNT);
        final SimpleInterval interval = new SimpleInterval(contig, start, end);
        return new SimpleCount(interval, count);
    };

    private static final BiConsumer<SimpleCount, DataLine> SIMPLE_COUNT_RECORD_TO_DATA_LINE_ENCODER = (simpleCount, dataLine) ->
            dataLine.append(simpleCount.getInterval().getContig())
                    .append(simpleCount.getInterval().getStart())
                    .append(simpleCount.getInterval().getEnd())
                    .append(simpleCount.getCount());

    private SimpleCountCollection(final File inputFile) {
        super(inputFile, SimpleCountCollection.SimpleCountTableColumn.COLUMNS, SIMPLE_COUNT_RECORD_FROM_DATA_LINE_DECODER, SIMPLE_COUNT_RECORD_TO_DATA_LINE_ENCODER);
    }

    public SimpleCountCollection(final SampleLocatableMetadata metadata,
                                 final List<SimpleCount> simpleCounts) {
        super(metadata, simpleCounts, SimpleCountCollection.SimpleCountTableColumn.COLUMNS, SIMPLE_COUNT_RECORD_FROM_DATA_LINE_DECODER, SIMPLE_COUNT_RECORD_TO_DATA_LINE_ENCODER);
    }

    /**
     * Read all counts from a file (HDF5 or TSV).
     */
    public static SimpleCountCollection read(final File file) {
        IOUtils.canReadFile(file);
        return readSubset(file, null);  //specifying a null intervalSubset returns all counts
    }

    /**
     * From a file (HDF5 or TSV), subset only the counts with intervals coinciding with intervals from a given list.
     * The list may contain intervals that do not coincide with any count intervals.
     * Unlike {@link #readSubset(String, List)}, this method first reads and constructs a {@link SimpleCountCollection}
     * using the entire file, and then creates and returns a second {@link SimpleCountCollection} containing only the
     * requested subset.
     * @param intervalSubset    if {@code null} or empty, all counts will be returned
     */
    public static SimpleCountCollection readSubset(final File file,
                                                   final List<SimpleInterval> intervalSubset) {
        IOUtils.canReadFile(file);
        final SimpleCountCollection simpleCounts = IOUtils.isHDF5File(file.toPath())
                ? readHDF5(new HDF5File(file))
                : readTSV(file);
        if (intervalSubset == null || intervalSubset.isEmpty()) {
            return simpleCounts;
        }
        return new SimpleCountCollection(
                simpleCounts.getMetadata(),
                simpleCounts.getRecords().stream()
                        .filter(c -> intervalSubset.contains(c.getInterval()))
                        .collect(Collectors.toList()));
    }

    private static SimpleCountCollection readTSV(final File file) {
        IOUtils.canReadFile(file);
        return new SimpleCountCollection(file);
    }

    private static SimpleCountCollection readHDF5(final HDF5File file) {
        Utils.nonNull(file);
        final HDF5SimpleCountCollection hdf5CountCollection = new HDF5SimpleCountCollection(file);
        final SampleLocatableMetadata metadata = hdf5CountCollection.getMetadata();
        final List<SimpleInterval> intervals = hdf5CountCollection.getIntervals();
        final double[] counts = hdf5CountCollection.getCounts().getRow(0);
        final List<SimpleCount> simpleCounts = IntStream.range(0, intervals.size())
                .mapToObj(i -> new SimpleCount(intervals.get(i), (int) counts[i]))
                .collect(Collectors.toList());
        return new SimpleCountCollection(metadata, simpleCounts);
    }

    /**
     * Read all counts from a Google Cloud Storage URL.
     */
    public static SimpleCountCollection read(final String path) {
        IOUtils.assertFileIsReadable(IOUtils.getPath(path));
        Utils.validate(BucketUtils.isCloudStorageUrl(path), "Read-count path must be a Google Cloud Storage URL.");
        return readSubset(path, null);  //specifying a null intervalSubset returns all counts
    }

    /**
     * From a Google Cloud Storage URL, subset only the counts with intervals coinciding with intervals from a given list.
     * The list may contain intervals that do not coincide with any count intervals.
     * @param intervalSubset    if {@code null} or empty, all counts will be returned
     */
    public static SimpleCountCollection readSubset(final String path,
                                                   final List<SimpleInterval> intervalSubset) {
        IOUtils.assertFileIsReadable(IOUtils.getPath(path));
        Utils.validate(BucketUtils.isCloudStorageUrl(path), "Read-count path must be a Google Cloud Storage URL.");
        final SAMTextHeaderCodec samTextHeaderCodec = new SAMTextHeaderCodec();
        final FeatureDataSource<SimpleCount> simpleCountsFeatureDataSource = new FeatureDataSource<>(
                path,
                path,
                DEFAULT_FEATURE_QUERY_LOOKAHEAD_IN_BP,
                SimpleCount.class,
                ConfigFactory.getInstance().getGATKConfig().cloudPrefetchBuffer(),
                ConfigFactory.getInstance().getGATKConfig().cloudIndexPrefetchBuffer());
        simpleCountsFeatureDataSource.setIntervalsForTraversal(intervalSubset);
        final BufferedLineReader reader = new BufferedLineReader(BucketUtils.openFile(path));
        final SAMFileHeader header = samTextHeaderCodec.decode(reader, path);
        final SampleLocatableMetadata metadata = MetadataUtils.fromHeader(header, Metadata.Type.SAMPLE_LOCATABLE);
        final List<SimpleCount> simpleCounts = Lists.newArrayList(simpleCountsFeatureDataSource.iterator());
        return new SimpleCountCollection(metadata, simpleCounts);
    }

    public void writeHDF5(final File file) {
        Utils.nonNull(file);
        HDF5SimpleCountCollection.write(file, getMetadata(), getIntervals(), getCounts());
    }

    public double[] getCounts() {
        return getRecords().stream().mapToDouble(SimpleCount::getCount).toArray();
    }
}
