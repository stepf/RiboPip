module Ribopip
  # Storing all filenames in one place for easier maintenance, inside a hash:
  # { key1: ['suffix1', 'fileinfo1'], key2: ['suffix2', 'fileinfo2'], ...}
  module Filenames
    # Generic class
    class Filenames
      # base - file basename
      attr_reader :base

      # initializes empty string and hash
      def initialize
        @base = ''
        @namehash = {}
      end

      # gets filename
      #
      # key - hash key
      # i   - if set, templating will be used
      #
      # Returns string or nil
      def get(key, i = nil)
        fail "Filenames: Could not get '#{key}'" if @namehash[key.to_sym].nil?
        if i.nil? # standard behavior
          "#{@base}#{@namehash[key.to_sym][0]}"
        else # !i.nil?: use templating for error rate
          "#{@base}#{@namehash[key.to_sym][0] % i}"
        end
      end

      # sets filename
      #
      # key   - hash key
      # value - new filename
      #
      # Returns string or nil
      def set(key, value)
        @namehash[key.to_sym][0] = value
      end

      # get info
      #
      # key - hash key
      #
      # Returns string or nil
      def get_info(key)
        @namehash[key.to_sym][1]
      end
    end

    # Alignment pipeline; also defines bucketing-specific functions
    class Alignment < Filenames
      # init filenames for pipeline
      def initialize(filename)
        @basename = File.basename(filename, File.extname(filename))
        @dirname  = File.dirname(filename)
        @base     = "#{@dirname}/#{@basename}"
        tophat_dir = '.vs_genome'
        @namehash = {
          counts: [
            '.counts',
            '#Reads after each alignment step, inluding buckets'
          ],
          reads: ['.fastq', 'Input reads'],
          readsstat: ['.fastq.len', 'Length distribution of input reads'],
          filter: ['.prep.filter.fastq', 'Reads after filtering'],
          clip: ['.prep.clip.fastq', 'Reads after adapter clipping'],
          cliplog: ['.prep.clip.log', 'Cutadapt clipping log'],
          trim: ['.prep.trim.fastq', 'Reads after 5\' trimming'],
          ncrna: ['.ncrna.sam', 'Reads aligning to ncRNA reference'],
          ncrnadist: [
            '.ncrna.dist.sam',
            'Reads distribution along ncRNA reference'
          ],
          fp: [
            '.fp.fastq',
            'Footprints: Reads not aligning to ncRNA reference'
          ],
          fpstat: ['.fp.lens', 'Footprint length distribution'],
          fpstatplot: ['.fp.lens.pdf', 'Plotted Footprint length distribution'],
          buckets: ['.fp.buckets', 'Length distribution of buckets'],
          topout: ["#{tophat_dir}", 'Output dir for Tophat alignment'],
          unmapped: [
            "#{tophat_dir}/unmapped.bam",
            'Reads not aligning to genomic reference'
          ],
          mapped_all: [
            "#{tophat_dir}/accepted_hits.bam",
            'Reads aligning to genomic reference'
          ],
          mapped_all_star: [
            "#{tophat_dir}/accepted_hits.bamAligned.out.sam",
            'Reads aligning to genomic reference'
          ],
          mapped_merged: ['.vs_genome.bam', 'Merged mapping reads'],
          unmapped_merged: [
            '.vs_genome.unmapped.bam',
            'Merged non-mapping reads'
          ],
          mapped_uniq: [
            '.vs_genome.uni.bam',
            'Reads aligning to genomic reference exactly once'
          ],
          mapped_uniqsort: [
            '.vs_genome.uni.sort.bam',
            'Sorted reads aligning to genomic ref exactly once'
          ],
          mapped_uni: [
            '.vs_genome.uni.sort.%ierr.bam',
            'Reads aligning to genomic reference with %i errors'
          ],
          mapped_igvtdf: [
            '.vs_genome.uni.sort.tdf',
            'Precomputed tdf track for faster display in IGV'
          ],
          mapped_igvwig: [
            '.vs_genome.uni.sort.wig',
            'Precomputed tdf track for faster display in IGV'
          ],
          ft: ['.vs_genome.uni.ft.txt', 'FeatureCounts Exon'],
          ftcds: ['.vs_genome.uni.ft.cds.txt', 'FeatureCounts CDS'],
          ftgenericpfx: ['.vs_genome.uni.ft.', 'FeatureCounts Generic Pfx'],
          ftgenericsfx: ['.txt', 'FeatureCounts Generic Sfx'],
          ftdist: [
            '.vs_genome.uni.ft.dist.txt',
            'FeatureCounts distribution'
          ],
          ftsums: [
            '.vs_genome.uni.ft.dist.sums.txt',
            'FeatureCounts distribution'
          ],
          ftnames: [
            '.vs_genome.uni.ft.dist.names.txt',
            'Feature Counts distribution with gene names'
          ],
          ftlog: [
            '.vs_genome.uni.ft.txt.summary',
            'Feature Counts summary'
          ],
          ftcdslog: [
            '.vs_genome.uni.ft.cds.txt.summary',
            'FeatureCounts CDS Log'
          ],
          ftgenericlogpfx: [
            '.vs_genome.uni.ft.',
            'FeatureCounts Generic Log Pfx'
          ]
          ftgenericlogsfx: [
            '.txt.summary',
            'FeatureCounts Generic Log Sfx'
          ]
        }
      end

      # Updates internal filenames according to current bucket
      #
      # lower - lower bound to read length
      # upper - upper bound to read length
      #
      # Returns nothing
      def set_bucket(lower, upper)
        @base = "#{@dirname}/#{lower}.#{upper}.#{@basename}"
      end

      # resets @base
      def unset_bucket
        @base = "#{@dirname}/#{@basename}"
      end
    end

    # Postprocessing Pipeline: DEseq
    class Deseq < Filenames
      # treatment - treatment name
      # ctrl      - control name
      attr_reader :treatment, :ctrl

      # init filenames for postprocessing
      #
      # treatment - treatment name
      # ctrl      - control name
      def initialize(treatment, ctrl, outdir)
        @base = "#{outdir}/#{treatment}_vs_#{ctrl}"
        @namehash = {
          merged: [
            '.deseq.counts.tsv', 'Merged FeatureCounts of all replicates'
          ],
          mergednorm: [
            '.deseq.counts.norm.tsv',
            'Merged and normalized FeatureCounts of all replicates'
          ],
          all: ['.deseq%s.all.csv', 'All DESeq%s results'],
          padj: ['.deseq%s.padj.csv', 'DESeq%s results with padj > 0.05'],
          log: ['.deseq%s.log', 'DESeq%s log'],
          allft: [
            '.deseq%s.all.ft.csv',
            'All DESeq%s results including raw FeatureCounts'
          ],
          padjft: [
            '.deseq%s.padj.ft.csv',
            'DESeq%s results with padj > 0.05 including raw FeatureCounts'
          ],
          allftnorm: [
            '.deseq%s.all.ftn.csv',
            'All DESeq%s results including raw FeatureCounts'
          ],
          padjftnorm: [
            '.deseq%s.padj.ftn.csv',
            'DESeq%s results with padj > 0.05 including raw FeatureCounts'
          ],
          allnames: [
            '.deseq%s.all.ft.names.csv',
            'All DESeq%s results including raw FeatureCounts'
          ],
          padjnames: [
            '.deseq%s.padj.ft.names.csv',
            'DESeq%s results with padj > 0.05 including raw FeatureCounts'
          ],
          xml: [
            '.deseq%s.padj.ft.comp.names.xml',
            'DESeq%s results with gene names'
          ]
        }
      end
    end

    # Postprocessing Pipeline: Further DESeq Processing
    class Postdeseq < Filenames
      # init filenames for postprocessing
      #
      # basename - input file
      def initialize(basename)
        @base = basename
        @namehash = {
          allft: [
            '.deseq%i.all.ft.csv',
            'All DESeq%i results including raw FeatureCounts'
          ],
          allftnorm: [
            '.deseq%i.all.ftn.csv',
            'All DESeq%i results including raw FeatureCounts'
          ],
          padjft: [
            '.deseq%i.padj.ft.csv',
            'DESeq%i results with padj > 0.05 including raw FeatureCounts'
          ],
          padjftnorm: [
            '.deseq%i.padj.ftn.csv',
            'DESeq%i results with padj > 0.05 including raw FeatureCounts'
          ],
          padjcp: [
            '.deseq%i.padj.ft.comp.csv',
            'DESeq%i results with padj > 0.05 including raw FeatureCounts'
          ],
          padjnames: [
            '.deseq%i.padj.ft.comp.names.csv',
            'DESeq%i results with gene names'
          ],
          padjxml: [
            '.deseq%i.padj.ft.comp.names.xml',
            'DESeq%i results with gene names'
          ]
        }
      end
    end

    # Postprocessing Pipeline: RiboDiff
    class Ribodiff < Filenames
      # init filenames for postprocessing
      #
      # treatment - treatment name
      # ctrl      - control name
      def initialize(treatment, ctrl, outdir)
        @base = "#{outdir}/#{treatment}_vs_#{ctrl}"
        @namehash = {
          merged: [
            '.rdiff.counts.tsv', 'Merged FeaterCounts of all replicates'
          ],
          outline: ['.rdiff.outline.csv', 'Experimental outline'],
          results: ['.rdiff.results.tsv', 'Raw results'],
          all: ['.rdiff.results.all.tsv', 'padj-filtered results'],
          padj: ['.rdiff.results.padj.tsv', 'padj-filtered results'],
          padjft: ['.rdiff.results.padj.ft.tsv', 'padj-filtered results'],
          padjnames: [
            '.rdiff.results.padj.ft.names.tsv',
            'padj-filtered results'
          ],
          xml: ['.rdiff.results.padj.ft.names.xml', 'padj-filtered results'],
          log: ['.rdiff.log', 'RiboDiff log']
        }
      end
    end

    # Postprocessing Pipeline: DESEq with TE
    class TE < Filenames
      # initializes names for translational efficiency pipeline
      def initialize(base)
        @base = "#{base}.TE.".gsub('fp_', '')
        @namehash = {
          counts: ['counts.tsv', 'Countains TE counts'],
          all: ['deseq1.all.tsv', 'Contains all results'],
          padj: ['deseq1.padj.tsv', 'Contains padj-filtered results'],
          padjft: ['deseq1.padj.ft.tsv', 'Contains padj-filtered results'],
          padjftfp: ['deseq1.padj.ftf.tsv', 'Contains padj-filtered results'],
          padjftmrna: [
            'deseq1.padj.ftfm.tsv',
            'Contains padj-filtered results'
          ],
          padjnames: [
            'deseq1.padj.ft.names.tsv',
            'Contains padj-filtered results'
          ],
          xml: ['deseq1.padj.ft.names.xml', 'Contains padj-filtered results'],
          log: ['deseq1.log', 'Contains DESeq output']
        }
      end
    end
  end
end
