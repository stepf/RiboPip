require 'ribopip'

module Ribopip
  # Implementing several strategies for differential expression analysis
  module ExpressionAnalysis
    # Making .write_rep_header available to ExpressionAnalysisStep
    # and AddTotalCountsToResults
    class RepHeaderGeneric
      # Generic function to compute and write header for merged featurecounts
      # according to condition names and numbers
      #
      # outfile    - File object
      # reps       - array of array strings; contains filename of each replicate
      #              e.g. [['file1', 'file2'], ['file3', 'file4']]
      # conditions - array of strings; contains names of experimental conditions
      def write_rep_header(outfile, reps, conditions)
        rep_names = []
        reps.each_with_index do |rep, i|
          1.upto(rep.length) { |j| rep_names << "#{conditions[i]}_rep#{j}" }
        end

        outfile.printf('%s', 'Entry')
        outfile.printf rep_names.map { "\t%s" }.join, *rep_names
        outfile.printf("\n")
      end
    end

    # Generic class for Deseq / Ribodiff Analysis
    class ExpressionAnalysisStep < RepHeaderGeneric
      include Ribopip
      # names      - object of class Filename; contains filenames & folders
      # features   - hash; key-value pairs of gene IDs to counts
      # reps       - array; filenames of replicates
      # reps_num   - integer; number of replicates
      # conditions - array; names of experimental conditions
      attr_reader :names, :features, :reps, :reps_num, :conditions

      # Initializes filenames, replicates, condition names and output folder
      #
      # names      - object of class Filename; contains filenames & folders
      # reps       - array of array strings; contains filename of each replicate
      #              e.g. [['file1', 'file2'], ['file3', 'file4']]
      # conditions - array of strings; contains names of experimental conditions
      # outdir     - output directory
      def initialize(names, reps, conditions, outdir)
        @names      = names
        @reps       = reps
        @conditions = conditions
        @outdir     = outdir
        @features   = {}
        @reps_num   = []
        @reps.each { |rep| @reps_num << rep.length }
      end

      # Wrapper function to merge counts of all replicates
      def merge_featurecounts
        replicates = @reps.flatten
        replicates.each_with_index do |rep, column|
          read_featurecounts(rep, column)
        end
        write_merged_featurecounts
      end

      # Generic function to parse featurecounts from tab-delimited file
      # into @features hash
      #
      # file   - input file
      # column - target colum in
      def read_featurecounts(file, column)
        File.readlines(file).each do |line|
          key = line.split[0]
          value = line.split[2].to_i
          @features[key.to_sym] = [] if @features[key.to_sym].nil? # init array
          # cumulative counting
          if @features[key.to_sym][column].nil?
            @features[key.to_sym][column] = value
          else
            @features[key.to_sym][column] += value
          end
        end
      end

      # Writes @features into tab-delimited file
      #
      # file - output file
      def write_merged_featurecounts
        file_out = File.open(@names.get('merged'), 'w')

        write_rep_header(file_out, @reps, @conditions)
        total_len = @reps.flatten.length
        @features.each do |gene_id, counts|
          file_out.printf("%s\t%s\n", gene_id, counts * "\t") \
            if counts.compact.length == total_len
        end
        file_out.close
      end
    end

    # Running RiboDiff
    class Ribodiff < ExpressionAnalysisStep
      # Runs ribodiff
      def compute
        compute_exp_outline
        merge_featurecounts
        run_cmd('ribodiff' \
                ' -d 0 -r 1 -p 1 -q 0.1' \
                " -e #{@names.get('outline')}" \
                " -c #{@names.get('merged')}" \
                " -o #{@names.get('results')}" \
                " > #{@names.get('log')}")
        filter_results
      end

      # Computes experimental outline from @reps and @conditions
      def compute_exp_outline
        rep_names = []
        cond_types = []
        cond_names = []
        @reps.each_with_index do |rep, i|
          1.upto(rep.length) do |j|
            rep_names << "#{@conditions[i]}_rep#{j}"
            # first 2 conditions are 'Ribo-Seq'
            cond_types << (i < 2 ? 'Ribo-Seq' : 'RNA-Seq')
            # every 2nd condition is 'Control'
            cond_names << (i.even? ? 'Treated' : 'Control')
          end
        end

        fail if (rep_names.length != cond_types.length) ||
                (rep_names.length != cond_names.length)

        file_out = File.open(@names.get('outline'), 'w')
        file_out.puts('Samples,Data_Type,Conditions')
        rep_names.each_with_index do |rep, i|
          file_out.printf("%s,%s,%s\n", rep, cond_types[i], cond_names[i])
        end
        file_out.close
      end

      # Filters @names.get('results'); writes entries with padj <= 0.05 into
      # @names.get('padj')
      def filter_results
        file_out = File.open(@names.get('padj'), 'w')

        # write header
        file_out.puts File.open(@names.get('results'), &:readline)
        # write body
        File.readlines(@names.get('results')).drop(1).each do |line|
          if line.split[3].to_f <= 0.05 && line.split[3] != 'nan'
            file_out.puts line
          end
        end

        file_out.close
      end
    end
  end
end
