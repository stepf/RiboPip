require 'ribopip'
require 'ribopip/array_writer'
require 'ribopip/expression_analysis'

module Ribopip
  # Implementing several strategies for differential expression analysis
  module ExpressionAnalysis
    # Using DESeq 1/2 for differential expression analysis
    module Deseq
      # Runs DESeq computation
      class Deseq < ExpressionAnalysisStep
        # Runs DESeq1 / 2 computation
        def compute
          merge_featurecounts
          run_deseq(1)
          run_deseq(2)
        end

        # Runs deseq
        #
        # ver - deseq version
        def run_deseq(ver)
          deseq_path = File.expand_path(
            '../../../scripts/run_deseq.R', __FILE__
          )
          cmd = "#{deseq_path} #{ver} #{@names.get('merged')} 128 " \
                  "#{@conditions[0]} #{@conditions[1]} #{@reps_num.join(' ')} "\
                  "#{@names.get('all', ver)} #{@names.get('padj', ver)} " \
                  "> #{@names.get('log', ver)}"
          run_cmd(cmd)
        end
      end

      # Identifies significant genes along both fp and total and includes them
      # in every other DESeq / DESeq2 report for comparison
      class ImportantGenes
        # Initializes arrays to contain IDs of significant genes
        #
        # basenames - basenames of deseq results
        def initialize(basenames)
          @basenames = basenames
          @significant_genes_deseq = []
          @significant_genes_deseq[1] = []
          @significant_genes_deseq[2] = []
          @all_genes_deseq = []
          @all_genes_deseq[1] = []
          @all_genes_deseq[2] = []
        end

        # Performs identification of significant genes, parses signifcant gene
        # IDS and writes them to output file
        def compute
          identify
          parse_all
          write
        end

        # Parses significant gene IDs from files in @basenames, stores them
        # separately in an corresponding array
        def identify
          @basenames.each_with_index do |names, idx|
            @significant_genes_deseq[1][idx] = []
            @significant_genes_deseq[2][idx] = []

            File.readlines("#{names.get('padjftnorm', 1)}").drop(1).each do |line|
              @significant_genes_deseq[1][idx].push(line.split[1].delete('"')) \
                unless line.split[1].nil?
            end

            File.readlines("#{names.get('padjftnorm', 2)}").drop(1).each do |line|
              @significant_genes_deseq[2][idx].push(line.split[0].delete('"'))\
                unless line.split[0].nil?
            end
          end
        end

        # Parses results for all genes from files in @basenames, stores them in
        # corresponding hashes for FP and total replicates
        def parse_all
          @basenames.each_with_index do |names, idx|
            @all_genes_deseq[1][idx] = {}
            @all_genes_deseq[2][idx] = {}

            (1..2).each do |ver|
              File.readlines(
                "#{names.get('allftnorm', ver)}"
              ).drop(ver).each do |line|
                key = line.split[ver].delete('"')
                @all_genes_deseq[ver][idx][key.to_sym] = line
              end
            end
          end
        end

        # Adds non-significant genes to results if they were significant in the
        # other experiment
        #
        # manual - manual targets
        def write(manual = [])
          @basenames.each_with_index do |names, idx|
            # Compute symmetric difference
            idx2 = idx == 0 ? 1 : 0
            targets_deseq = []
            targets_deseq[1] =
              @significant_genes_deseq[1][idx2].to_a - \
              @significant_genes_deseq[1][idx].to_a
            targets_deseq[2] =
              @significant_genes_deseq[2][idx2].to_a - \
              @significant_genes_deseq[2][idx].to_a

            # Write to files
            (1..2).each do |ver|
              File.open("#{names.get('padjcp', ver)}", 'a') do |file|
                file.puts # newline
                (targets_deseq[ver] + manual).each do |target|
                  file.puts @all_genes_deseq[ver][idx][target.to_sym] \
                    unless @all_genes_deseq[ver][idx][target.to_sym].nil?
                end
              end
            end
          end
        end
      end

      # Parses tab-delimited counts table into hash and adds them to any
      # DESeq table
      class AddCountsToResults
        # hash containing gene IDs as keys and counts as values
        attr_reader :features
        # Initialize infile and suffix; parses counts from infile into hash
        #
        # infile       - input file
        # suffix       - suffix to add to replicate names
        # norm_factors - DESeq log; if != nil normalized counts will be written
        def initialize(infile, suffix = nil, deseq_log = nil)
          @features = parse_counts(infile)
          @rep_names = parse_rep_names(infile)
          @rep_names.length.times { |i| @rep_names[i] += suffix } \
            unless suffix.nil?
          @norm_factors = deseq_log.nil? ? nil : parse_norm_factors(deseq_log)
        end

        # Reads DESeq results and adds parsed counts to them
        #
        # column  - col number containing the key (gene id)
        # infile  - input file
        # outfile - output file
        def compute(column, infile, outfile)
          header = File.open(infile, &:readline).split + @rep_names
          file_out = File.open(outfile, 'w')
          file_out.printf(header.map { "\t%s" }.join, *header)
          file_out.printf("\n")
          File.readlines(infile).drop(1).each do |line|
            key = line.split[column].delete('"')
            features = @norm_factors.nil? ? @features[key.to_sym] : normalize_counts(@features[key.to_sym])
            line_new = line.split + features
            file_out.puts line_new.flatten.join("\t")
          end
          file_out.close
        end

        # Writes normalized counts into new file for later TE computation
        #
        # outfile - output file
        def write(outfile)
          out = File.open(outfile, 'w')
          # write header
          header = @rep_names
          out.printf(header.map { "\t%s" }.join, *header)
          out.printf("\n")
          # write body
          @features.each do |key, counts|
            norm_counts = normalize_counts(counts)
            out.puts [key, norm_counts].flatten.join("\t")
          end
          out.close
        end

        # Parses counts into hash
        #
        # infile - input file
        #
        # Returns hash
        def parse_counts(infile)
          features = {}
          File.readlines(infile).drop(1).each do |line|
            key = line.split[0]
            values = line.split[1..-1]
            features[key.to_sym] = values.map(&:to_i)
          end
          features
        end

        # Parses replicate names from input file
        #
        # infile - input file
        #
        # Returns array of strings
        def parse_rep_names(infile)
          File.open(infile, &:readline).split[1..-1]
        end

        # Parses normalization factors from Deseq1 log file
        #
        # infile - input file
        #
        # Returns array of integers
        def parse_norm_factors(infile)
          line = 5 # line number that contains the size factors
          norm_factors = []
          # until total number of replicates is reached
          until norm_factors.length == @rep_names.length
            norm_factors << File.readlines(infile)[line].split.map(&:to_f)
            norm_factors.flatten! # otherwise .length estimation will fail
            line += 2 # in case norm factors span multiple lines
          end
          norm_factors.flatten
        end

        # Normalizes counts according to normalization factors in norm_factors
        #
        # counts - input array of integers
        #
        # Returns array of integers
        def normalize_counts(counts)
          normalized = []
          counts.each_with_index do |count, index|
            normalized << count * @norm_factors[index]
          end
          normalized
        end
      end

      # parses (normalized) counts from FP and mRNA and computes translational
      # efficiency by dividing counts(fp) / counts (mrna)
      class TE
        include Ribopip

        # init parameters; TODO: seperate .initialize method
        #
        # fp      - path to fp counts
        # mrna    - path to mrna counts
        # outfile - output file
        def self.compute(fp, mrna, outfile)
          @fp_counts = read_counts(fp)
          @mrna_counts = read_counts(mrna)
          @te = {}
          compute_te
          write_te(outfile, fp)
        end

        # read counts from tab-delimited file
        #
        # infile - input file
        def self.read_counts(infile)
          counts = {}
          File.readlines(infile).drop(1).each do |line|
            key = line.split[0]
            values = line.split[1..-1].map(&:to_f)
            counts[key.to_sym] = values
          end
          counts
        end

        # computs translational efficiency using a workaround:
        #   TE values are multiplied by 1000 and then rounded as DESeq does not
        #   support floats as input values
        def self.compute_te
          @fp_counts.each do |key, values|
            unless @mrna_counts[key.to_sym].nil?
              te_values = []
              values.each_with_index do |value, i|
                te_values << (value * 1000 / @mrna_counts[key.to_sym][i]).round
              end
            end
            @te[key.to_sym] = te_values
          end
        end

        # writes TE values to outfile
        #
        # outfile - output file
        # fp      - copy header from one of the input files
        def self.write_te(outfile, fp)
          file_out = File.open(outfile, 'w')
          file_out.puts File.open(fp, &:readline) # write header
          num = 0
          @te.each do |gene_id, counts| # write genes
            if !counts.nil? && !counts.empty?
              file_out.puts [gene_id, counts].flatten.join("\t")
            else
              num += 1
            end
          end
          file_out.close
        end
      end

      # Converts Deseq results to Excel-compatible XML
      class Deseq2Excel
        # read file into array, then convert
        #
        # infile  - input file
        # outfile - output file
        def self.convert(infile, outfile)
          xml_array = []

          File.readlines(infile).each do |line|
            xml_array << line.chomp.split(/\t/)
          end

          writer = Ribopip::ArrayWriter::XMLWriter.new(xml_array)
          writer.write(outfile)
        end
      end
    end
  end
end
