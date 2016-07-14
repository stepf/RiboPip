require 'ribopip'

module Ribopip
  # Pre-Processing and Alignment Pipeline for Ribosome Profiling Data.
  module Pipeline
    # abstract class for a pipeline step
    class PipelineStep
      include Ribopip

      # Initializes all output filenames and folders for later use
      #
      # names           - object of class Filename; contains filenames & folders
      # force_overwrite - If true, pipeline will overwrite all existing files
      #
      # Returns nothing
      def initialize(names, force_overwrite)
        @names           = names
        @force_overwrite = force_overwrite
      end

      # Verbosely checks if file should be processed
      #
      # filename - file whose existance will be checked
      # step     - name of pipeline step
      #
      # Returns true if file does not exist or @force_overwrite is set true.
      # Otherwise returns false.
      def skip_step?(filename, step)
        if File.exist?(filename) && !@force_overwrite
          print_e "SKIPPED #{step}: #{filename} already exists."
          true
        else
          print_e "RUN #{step} => #{filename}"
          false
        end
      end
    end

    # Pre-Processing steps before actual alignment
    class Preprocessing < PipelineStep
      # Initializes Prep-Parameters
      #
      # minlen   - discard all shorter reads
      # linker   - string to be clipped from 3' ends of reads
      # software - clipping software
      #            :fastx    - FASTX Toolkit
      #            :cutadapt - Cutadapt
      # err_rate - allowed error rate for clipping (only cutadapt)
      def initialize(names, force_overwrite, minlen, linker, software, err_rate)
        super(names, force_overwrite)
        @minlen = minlen
        @linker = linker
        @software = software
        @err_rate = err_rate
      end

      # Performs Pre-Processing
      def compute
        filter(@minlen)
        clip(@linker, @software, @err_rate, @minlen)
        trim
      end

      # Quality filters reads
      #
      # minlen - discard all shorter reads
      #
      # Returns nothing
      def filter(minlen)
        return if skip_step?(@names.get('filter'), 'filtering')

        # Only filter input files from Illumina CASAVA 1.8 pipeline
        if `head -n 1 #{@names.get('reads')} | cut -d ' ' -f 3`.empty?
          run_cmd(
            'fastq_illumina_filter' \
            " --keep N -v -l #{minlen} " \
            " -o #{@names.get('filter')}" \
            " #{@names.get('reads')}"
          )
        else
          @names.set('filter', '.fastq')
        end
      end

      # Clips linker from 3'
      #
      # linker     - string to be clipped from 3'
      # software   - clipping software (fastx, cutadapt)
      # error_rate - allowed error rate (only for cutadapt)
      # minlen     - discard all shorter reads
      #
      # Returns nothing.
      def clip(linker, software, error_rate, minlen)
        return if skip_step?(@names.get('clip'), 'clipping')

        clipper_cmd = {
          fastx: \
            'fastx_clipper' \
            ' -Q33 -c -n -v' \
            " -a #{linker}" \
            " -l #{minlen}" \
            " #{@names.base}" \
            " -i #{@names.get('filter')}" \
            " -o #{@names.get('clip')}",
          cutadapt: \
            'cutadapt' \
            " -a #{linker}" \
            ' --trimmed-only' \
            " -e #{error_rate}" \
            " -m #{minlen}" \
            " -o #{@names.get('clip')}" \
            " #{@names.get('filter')}" \
            "> #{@names.get('cliplog')}"
        }
        run_cmd(clipper_cmd[software])
      end

      # Trims nucleotide from the 5' end of each read
      #
      # Returns nothing
      def trim
        return if skip_step?(@names.get('trim'), 'trimming')
        run_cmd(\
          'fastx_trimmer -Q33 -f 2' \
          " -i #{@names.get('clip')}" \
          " -o #{@names.get('trim')}"
        )
      end
    end # class Preprocessing

    # Abstract class for alignment step; both ncRNA and genomic
    class AlignmentPipelineStep < PipelineStep
      # Initializes alignment parameters
      #
      # ref        - path to reference
      # software   - alignment software (bowtie1, bowtie2, bwa, star)
      def initialize(names, force_overwrite, ref, software)
        super(names, force_overwrite)
        @ref = ref
        @software = software
        @ref_base = "#{File.dirname(ref)}/#{File.basename(ref, '.fa')}"
      end

      # Pre-Computes index if it does not exist
      #
      # ref        - path to reference
      # ref_base   - path to reference without file extension
      # software   - alignment software (bowtie1, bowtie2, bwa, star)
      # annotation - path to GTF annotation (only star)
      #
      # Returns nothing
      def index(ref, ref_base, software, annotation = '')
        index_suffix = {
          bowtie1: '4.ebwt',
          bowtie2: '4.bt2',
          bwa: '.sa',
          star: '.star'
        }
        index_cmd = {
          bowtie1: "bowtie-build -p #{ref} #{ref_base}",
          bowtie2: "bowtie2-build -p #{ref} #{ref_base}",
          bwa: "bwa index #{ref}",
          star: "mkdir #{ref_base} && "\
            'STAR --runMode genomeGenerate' \
            ' --runThreadN $(nproc)' \
            " --genomeDir #{ref_base}"\
            " --genomeFastaFiles #{ref}"\
            ' --sjdbOverhang 49' \
            " --sjdbGTFfile #{annotation} "
        }

        time = 5
        while File.exist?("#{ref_base}.lock")
          print_e "#{ref_base}.lock exists. Wait for #{time} seconds."
          sleep(time)
          time *= 5
        end

        return if software == :tophat ||
                  skip_step?("#{ref_base}.#{index_suffix[software]}", 'indexing')

        begin
          run_cmd("touch #{ref_base}.lock")
          run_cmd(index_cmd[software])
        ensure
          run_cmd("rm -f #{ref_base}.lock")
        end
      end

      # Performs alignment
      #
      # ref      - path to genomic reference
      # ref_base - path to reference without file extension
      # software - alignment software (bowtie1, bowtie2, bwa, star, tophat)
      # opts     - step-specific options
      #            :annotation     - path to genomic annotation
      #            :mismatches     - max num of mismatches in alignment
      #            :seedlen        - seed length for ncRNA alignment
      #            :tophat_aligner - software that tophat uses (bowtiet1 / 2)
      #
      # Returns name of files containing mapped and unmapped reads
      def align(ref, ref_base, software, opts = {})
        if software == :tophat
          bt_flag =
            opts[:tophat_aligner] == :bowtie1 ? '--bowtie1' : ''
          gap_flag =
            opts[:mismatches] < 2 ? "--read-gap-length #{opts[:mismatches]}" : ''
        end

        aln_cmd = {
          bowtie1:
            'bowtie' \
            " --seedlen=#{opts[:seedlen]} #{ref_base}" \
            " --un=#{@names.get('fp')}" \
            " -q #{@names.get('trim')} " \
            " --sam #{@names.get('ncrna')}",
          bowtie2:
            'bowtie2' \
            " --un #{@names.get('fp')}" \
            " -x #{ref_base}" \
            " -L #{opts[:seedlen]}" \
            " -U #{@names.get('trim')}" \
            " -S #{@names.get('ncrna')}",
          bwa:
            'bwa mem' \
            " -k #{opts[:seedlen]}" \
            " #{ref} " \
            " #{@names.get('trim')} " \
            "| samtools view -b - > #{@names.get('ncrna')} " \
            '&& bam2fastq' \
            " -o #{@names.get('fp')}" \
            " --no-aligned #{@names.get('ncrna')}",
          tophat:
            'tophat' \
            " --read-edit-dist #{opts[:mismatches]}" \
            " #{bt_flag}" \
            " -N #{opts[:mismatches]}" \
            " --output-dir #{@names.get('topout')}" \
            ' --no-novel-juncs' \
            " #{gap_flag}" \
            " --GTF #{opts[:annotation]}" \
            " #{ref_base} #{@names.get('fp')}",
          star:
            'STAR' \
            " --genomeDir #{ref_base}" \
            " --outFilterMismatchNmax #{opts[:mismatches]}" \
            " --readFilesIn #{@names.get('fp')}"\
            " --outFileNamePrefix #{@names.get('mapped_all')}"
        }

        target =
          opts[:seedlen].nil? ? @names.get('mapped_all') : @names.get('fp')
        run_cmd(aln_cmd[software]) unless skip_step?(target, 'aligning')
        [@names.get('mapped_all'), @names.get('unmapped')]
      end
    end # class AlignmentPipelineStep

    # Alignment to ncRNA reference
    class NcrnaAlignment < AlignmentPipelineStep
      # Initializes ncRNA alignment parameters
      #
      # seedlen - seed length
      def initialize(names, force_overwrite, ref, software, seedlen)
        super(names, force_overwrite, ref, software)
        @seedlen = seedlen
      end

      # Performs alignment to ncRNA reference; only unaligned reads will be
      # processed further
      def compute
        index(@ref, @ref_base, @software)
        align(@ref, @ref_base, @software, {seedlen: @seedlen})
      end
    end # class NcrnaAlignment

    # Splice-aware alignment to genomic reference and annotation
    class GenomicAlignment < AlignmentPipelineStep
      # max_mismatches - integer; number of allowed mismatches
      attr_reader :max_mismatches

      # Initializes genomic alignment parameters
      #
      # annotation     - genomic annotation, needed for splicing awareness
      # tophat_aligner - software that Tophat uses (bowtie2 or bowtie1)
      # mismatches     - max number of allowed mismatches in Tophat alignment
      # err            - error rate for STAR / bucketizing
      def initialize(names, force_overwrite, ref, software,
                     annotation, tophat_aligner, mismatches, err_rate)
        super(names, force_overwrite, ref, software)
        @annotation = annotation
        @tophat_aligner = tophat_aligner
        @mismatches = mismatches
        @err_rate = err_rate
        @mapped_bams = []
        @unmapped_bams = []
        @max_mismatches = 0
      end

      # Performs genomic alignment. As Tophat does only allow for an absolute
      # number of errors, reads will be split into several files according to
      # their length, if @err_rate is set. Then each file (bucket) will be
      # aligned seperately and the resulting alignments will be merged.
      def compute
        index(@ref, @ref_base, @software, @annotation)

        if @err_rate > 0 # TODO: evaluate && @software == :tophat
          bucketized_alignment
        else # software == :star || err_rate == 0
          unbucketized_alignment
        end
      end

      # Performs genomic alignment with relative error rate (buckets)
      def bucketized_alignment
        # split reads into buckets according to their size and err_rate
        @buckets = bucketize(@err_rate)

        # perform alignment on each bucket
        @buckets.reverse_each do |lower, upper, mismatches|
          @names.set_bucket(lower, upper)
          mapped, unmapped = align(
            @ref, @ref_base, @software,
            { annotation: @annotation,
              tophat_aligner: @tophat_aligner,
              mismatches: mismatches
            }
          )
          @mapped_bams << mapped
          @unmapped_bams << unmapped
          @max_mismatches = [@max_mismatches, mismatches].max
        end

        # merge alignments
        @names.unset_bucket
        unbucketize(@mapped_bams, @names.get('mapped_merged'))
        unbucketize(@unmapped_bams, @names.get('unmapped_merged'))
      end

      # Performs genomic alignment with fixed number of allowed mismatches
      def unbucketized_alignment
        align(
          @ref, @ref_base, @software,
          { annotation: @annotation,
            tophat_aligner: @tophat_aligner,
            mismatches: @mismatches
          }
        )
        mapped_all = @software == :star ? \
          @names.get('mapped_all_star') : @names.get('mapped_all')
        run_cmd("cp #{mapped_all} #{@names.get('mapped_merged')}")
        unless @software == :star
          run_cmd(
            "cp #{@names.get('unmapped')} #{@names.get('unmapped_merged')}"
          )
        end
        @max_mismatches = @mismatches
      end

      # Splits reads into several files containing a read length range to allow
      # for seperate alignments with a relative number of errors.
      #
      # error rate - = num of errors / read length
      #
      # Returns array.
      def bucketize(error_rate)
        buckets = []
        run_cmd(
          "fastq-bucketize #{@names.get('fp')} #{error_rate} " \
          "2> #{@names.get('buckets')}"
        )

        # parse buckets and compute corresponding absolute number of errors
        File.readlines(@names.get('buckets')).each do |line|
          next if line[0] == '#' # comment
          line = line.split.map(&:to_i)
          fail if (line[0] * error_rate).floor != (line[1] * error_rate).floor
          # push [lower bound, upper bound, absolute #errors]
          buckets.push([line[0], line[1], (line[0] * error_rate).floor]) \
            unless line[1] < 14 # TODO: implement minlen option
        end

        buckets
      end

      # Merges bucketed alignments into single bam file
      #
      # infiles - array of filenames
      # outfile - merged file
      def unbucketize(infiles, outfile)
        run_cmd("samtools merge -f #{outfile} #{infiles.join(' ')}")
      end
    end # class GenomicAlignment

    # Extract reads from alignment according to their #mismatches
    class ReadExtractor < PipelineStep
      # Initializes read extraction parameters
      #
      # mis - max number of mismatches
      def initialize(names, force_overwrite, mis)
        super(names, force_overwrite)
        @mis = mis
      end

      # Performs extraction
      def compute
        #return if skip_step?(@names.get('mapped_uni', @mis), 'read extraction')
        extract_uni
        (0..@mis).each do |i|
          extract_uni_err(i)
        end
      end

      # Extracts uniquely mapping reads
      def extract_uni
        # Extract all uniquely mapping reads
        run_cmd(
          "samtools view -H #{@names.get('mapped_merged')} " \
          "> #{@names.get('mapped_uniq')}"
        )
        run_cmd(
          "samtools view -h #{@names.get('mapped_merged')} " \
          "| grep -P 'NH:i:1(\t|$)' " \
          ">> #{@names.get('mapped_uniq')}"
        )
        run_cmd(
          "samtools sort -o #{@names.get('mapped_uniqsort')} -O bam -T " \
          "tmp.bam #{@names.get('mapped_uniq')}"
        )
      end

      # Extracts uniquely mapping reads with i mismatches
      #
      # mis - number of mismatches
      def extract_uni_err(mis)
        run_cmd(
          "samtools view -H #{@names.get('mapped_uniqsort')}" \
          "> #{@names.get('mapped_uni', mis)}"
        )
        run_cmd(
          "samtools view -h #{@names.get('mapped_uniqsort')} " \
          "| grep -E '([nN]M:i:#{mis})|(^@)' " \
          '| samtools view -S - ' \
          ">> #{@names.get('mapped_uni', mis)}"
        )
      end
    end # class ReadExtractor
  end # module Pipeline
end # module Ribopip
