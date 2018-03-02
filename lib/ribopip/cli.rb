require 'thor'
require 'ribopip'

# rubocop:disable all
module Ribopip
  # the command line interface
  class CLI < Thor
    include Ribopip

    desc 'align', 'runs alignment pipeline and generates corresponding metrics'
    method_option :reads,
      :aliases => '-r',
      :banner => 'FASTQ',
      :desc => 'Path to reads.',
      :required => true,
      :type => :string
    method_option :ncrna_ref,
      :aliases => '-n',
      :banner => 'FASTA',
      :desc => 'Path to ncRNA reference.',
      :required => true,
      :type => :string
    method_option :genomic_ref,
      :aliases => '-g',
      :banner => 'FASTA',
      :desc => 'Path to genomic reference.',
      :required => true,
      :type => :string
    method_option :genomic_annotation,
      :aliases => '-a',
      :banner => 'GTF',
      :desc => 'Path to genomic annotation.',
      :required => true,
      :type => :string
    method_option :clip_linker,
      :aliases => '-l',
      :banner => 'STRING',
      :default => 'CTGTAGGCACCATCAAT',
      :desc => 'String to be clipped from 3\'',
      :type => :string
    method_option :clip_software,
      :banner => 'PROGRAM',
      :enum => ['fastx', 'cutadapt'],
      :default => 'cutadapt',
      :desc => 'Clipping software',
      :type => :string
    method_option :clip_err_rate,
      :banner => 'FLOAT',
      :default => 0.2,
      :desc => 'Maximum allowed error rate (FLOAT = #errors / region length).',
      :type => :numeric
    method_option :clip_minlen,
      :banner => 'INT',
      :desc => 'Discard reads shorter than INT.',
      :default => 0,
      :type => :numeric
    method_option :ncrna_software,
      :default => 'bowtie1',
      :banner => 'PROGRAM',
      :desc => 'Software for ncRNA alignment',
      :enum => ['bowtie1', 'bowtie2', 'bwa'],
      :type => :string
    method_option :ncrna_seedlen,
      :banner => 'INT',
      :desc => 'Seed length for bowtie1 / bowtie2 alignment',
      :default => 14,
      :type => :numeric
    method_option :ncrna_minlen,
      :banner => 'INT',
      :desc => 'Discard reads shorter than INT',
      :default => 14,
      :type => :numeric
    method_option :genomic_software,
      :default => 'tophat',
      :banner => 'PROGRAM',
      :desc => 'Software for genomic alignment',
      :enum => ['tophat', 'star'],
      :type => :string
    method_option :genomic_tophat_engine,
      :default => 'bowtie1',
      :banner => 'PROGRAM',
      :desc => 'Software for genomic alignment using tophat',
      :enum => ['bowtie2', 'bowtie1'],
      :type => :string
    method_option :genomic_mismatches,
      :banner => 'INT',
      :desc => 'Max num of mismatches for genomic alignment',
      :default => 3,
      :type => :numeric
    method_option :genomic_err_rate,
      :banner => 'FLOAT',
      :desc => 'Set max num of mismatches for genomic alignment relative to ' \
        'read length (#errors = FLOAT * read length), if FLOAT > 0. ' \
        'Will overrule --genomic_mismatches.',
      :default => 0,
      :type => :numeric
    method_option :igv_ref,
      :banner => 'GENOME',
      :desc => 'Path to IGV .genome file. If set TDF track will be computed.',
      :type => :string
    method_option :force_overwrite,
      :desc => 'Overwrite all existing files',
      :default => false,
      :type => :boolean
    # runs alignment pipeline and generates corresponding metrics
    def align
      require 'ribopip/array_writer'
      require 'ribopip/binary_checker'
      require 'ribopip/counts'
      require 'ribopip/feature_counts'
      require 'ribopip/filenames'
      require 'ribopip/pipeline'
      require 'ribopip/metrics'

      bin_array = [
                    ['fastq_illumina_filter', '-h', '0.0.13', 'version'],
                    ['fastx_clipper', '-h', '0.0.13', 'version'],
                    ['fastx_trimmer', '-h', '0.0.13', 'version'],
                    ['cutadapt', '--version', '1.8.1', ''],
                    ['bowtie', '--version', '1.1.2', 'version'],
                    ['bowtie2', '--version', '2.2.5', 'version'],
                    ['bwa', '', '0.7.4', 'Version:'],
                    ['bam2fastq', '-v', '1.1.0', 'bam2fastq v'],
                    ['STAR', '', '2.4.1', 'STAR_'],
                    ['tophat', '--version', '2.0.13', 'TopHat v'],
                    ['fastq-bucketize', '-v', '.0.1', 'fastq-bucketize v.'],
                    ['samtools', '--version', '1.2', 'samtools'],
                    ['fastqc', '--version', '0.10.1', 'FastQC v'],
                    ['gt', '--version', '1.5.5', '(GenomeTools)'],
                    ['featureCounts', '-v', '1.4.4', 'featureCounts v'],
                    ['igvtools', 'version', '2.3.60', 'Version']
                  ]

      checker = Ribopip::BinaryChecker.new(bin_array)
      checker.perform_checks

      begin
        ##
        # INIT
        names = Ribopip::Filenames::Alignment.new(options[:reads])
        counts = Ribopip::Counts::Bookkeeper.new
        counts.parse_nreads(names.get('reads'), 'n_total')

        #Metrics::GtSeqstat.compute(names.get('reads'),
        #                              names.get('readsstat'))

        ##
        # PREPROCESSING: Quality filter, adapter clipping, 5' trimming
        prep = Ribopip::Pipeline::Preprocessing.new(
          names,
          options[:force_overwrite],
          options[:clip_minlen],
          options[:clip_linker],
          options[:clip_software].to_sym,
          options[:clip_err_rate]
        )
        prep.compute

        if options[:clipping_software] == 'cutadapt'
          counts.parse_nreads(
            names.get('cliplog'), 'n_postclip_adapter', 'with adapters'
          )
          counts.parse_nreads(
            names.get('cliplog'), 'n_postclip_length', 'passing filters'
          )
        end

        counts.parse_nreads(names.get('trim'), 'n_posttrim')

        ##
        # ALIGN TO NCRNA REFERENCE: Only unaligning reads are further processed
        ncrna = Ribopip::Pipeline::NcrnaAlignment.new(
          names,
          options[:force_overwrite],
          options[:ncrna_ref],
          options[:ncrna_software].to_sym,
          options[:ncrna_seedlen]
        )
        ncrna.compute
        counts.parse_nreads(names.get('ncrna'), 'n_ncrna_contaminants')
        counts.parse_nreads(names.get('fp'), 'n_reads_of_interest')

        # Compute ncRNA distribution over reference
        Metrics::FieldDistri.compute(
          names.get('ncrna'), names.get('ncrnadist'), options[:ncrna_ref]
        )
        # Compute and plot fp length distribution
        Metrics::GtSeqstat.compute(names.get('fp'), names.get('fpstat'))
        Metrics::GtSeqstat.plot(names.get('fpstat'), names.get('fpstatplot'))

        ##
        # ALIGN TO GENOMIC REF
        genomic =
          Ribopip::Pipeline::GenomicAlignment.new(
            names,
            options[:force_overwrite],
            options[:genomic_ref],
            options[:genomic_software].to_sym,
            options[:genomic_annotation],
            options[:genomic_tophat_engine].to_sym,
            options[:genomic_mismatches],
            options[:genomic_err_rate]
          )
        genomic.compute

        Metrics::FastQC.compute(names.get('mapped_merged'))
        Metrics::FastQC.compute(names.get('unmapped_merged'))
        counts.parse_nreads("#{names.get('mapped_merged')}", 'n_mapped')

        ##
        # EXTRACT UNIQUELY MAPPING READS WITH X MISMATCHES
        extraction =
          Ribopip::Pipeline::ReadExtractor.new(
            names, options[:force_overwrite], genomic.max_mismatches
          )
        extraction.compute

        counts.parse_nreads("#{names.get('mapped_uniq')}", 'n_mapped_uni')
        genomic.max_mismatches.downto(0) do |i|
          counts.parse_nreads(
            "#{names.get('mapped_uni', i)}", "n_mapped_uni_#{i}err"
          )
        end

        FeatureCounts.compute(
          "#{names.get('mapped_uniq')}",
          "#{names.get('ft')}",
          options[:genomic_annotation]
        )
        FeatureCounts.compute_distri(
          "#{names.get('ft')}",
          "#{names.get('ftdist')}",
          "#{names.get('ftsums')}",
          options[:genomic_annotation]
        )
        unless options[:igv_ref].nil?
          Metrics::IGV.compute(
            names.get('mapped_merged'),
            names.get('mapped_igvtdf'),
            names.get('mapped_igvwig'),
            options[:igv_ref]
          )
        end

        counts.parse_nreads("#{names.get('ftlog')}", 'n_assigned', 'Assigned')
        counts.parse_nreads(
          "#{names.get('ftsums')}", 'n_assigned_pc', 'protein_coding'
        )

        # get #reads mapping to CDS or UTR
        FeatureCounts.compute(
          "#{names.get('mapped_uniq')}",
          "#{names.get('ftcds')}",
          options[:genomic_annotation],
          'CDS'
        )
        FeatureCounts.compute(
          "#{names.get('mapped_uniq')}",
          "#{names.get('ftutr')}",
          options[:genomic_annotation],
          'UTR'
        )
        counts.parse_nreads("#{names.get('ftcdslog')}", 'n_cds', 'Assigned')
        counts.parse_nreads("#{names.get('ftutrlog')}", 'n_utr', 'Assigned')
      ensure
        ##
        # WRITE COUNTS
        tsv = Ribopip::ArrayWriter::TSVWriter.new(counts.array, true)
        tsv.write("#{names.get('counts')}.tsv")
        xml = Ribopip::ArrayWriter::XMLWriter.new(counts.array, true)
        xml.write("#{names.get('counts')}.xml")
      end
    end

    desc 'postproc', 'runs downstream analysis pipeline'
    method_option :annotation,
      :aliases => '-a',
      :banner => 'GTF',
      :desc => 'Path to genomic annotation',
      :required => true,
      :type => :string
    method_option :outdir,
      :aliases => '-o',
      :desc => 'Output directory',
      :required => true,
      :type => :string
    method_option :padj_threshold,
      :aliases => '-p',
      :banner => 'INT',
      :desc => 'Significance threshold for the adjusted p-value',
      :default => 0.05,
      :type => :numeric
    method_option :min_coverage,
      :aliases => '-m',
      :banner => 'INT',
      :desc => 'Minimum number of assignments to a gene across investigated replicates',
      :default => 128,
      :type => :numeric
    method_option :fp_1,
      :banner => 'NAME REP1 REP2 [REP3]',
      :desc => 'Treatment1 Replicates (e.g. mutant)',
      :type => :array
    method_option :fp_2,
      :banner => 'NAME REP1 REP2 [REP3]',
      :desc => 'Treatment2 Replicates (e.g. mutant)',
      :type => :array
    method_option :fp_control,
      :banner => 'NAME REP1 REP2 [REP3]',
      :desc => 'Control Replicates',
      :type => :array
    method_option :mrna_1,
      :banner => 'NAME REP1 REP2 [REP3]',
      :desc => 'Treatment1 Replicates (e.g. mutant)',
      :type => :array
    method_option :mrna_2,
      :banner => 'NAME REP1 REP2 [REP3]',
      :desc => 'Treatment2 Replicates (e.g. mutant)',
      :type => :array
    method_option :mrna_control,
      :banner => 'NAME REP1 REP2 [REP3]',
      :desc => 'Control Replicates',
      :type => :array
    method_option :ribodiff,
      :default => false,
      :desc => 'Run RiboDiff Analysis',
      :type => :boolean
    # runs downstream analysis pipeline
    def postproc
      require 'ribopip/deseq'
      require 'ribopip/expression_analysis'
      require 'ribopip/filenames'
      require 'ribopip/gene_id_2_name'

      ##
      # DESEQ
      fp_experiments =
        [options[:fp_1], options[:fp_2], options[:fp_control]].compact
      mrna_experiments  =
        [options[:mrna_1], options[:mrna_2], options[:mrna_control]].compact
      deseq_results = []

      # Run DESeq for each condition against each other condition
      [fp_experiments, mrna_experiments].each do |experiments|
        experiments.length.times do |i|
          (i+1).upto(experiments.length-1) do |j|
            condition1 = experiments[i]
            condition2 = experiments[j]
            names = Ribopip::Filenames::Deseq.new(
              condition1[0], condition2[0], options[:outdir]
            )
            deseq = Ribopip::ExpressionAnalysis::Deseq::Deseq.new(
              names.dup,
              [condition1[1..-1], condition2[1..-1]],
              [condition1[0],condition2[0]],
              options[:outdir]
            )
            deseq.compute
            deseq_results.push(deseq)
          end
        end
      end

      id2name = GeneID2Name.new(options[:annotation])
      deseq_results.each do |deseq|
        counts = Ribopip::ExpressionAnalysis::Deseq::AddCountsToResults.new(
          deseq.names.get('merged')
        )
        # add counts to unfiltered results
        counts.compute(
          1, deseq.names.get('all', 1), deseq.names.get('allft', 1)
        )
        counts.compute(
          0, deseq.names.get('all', 2), deseq.names.get('allft', 2)
        )
        # add counts to significant results (p-value)
        counts.compute(
          1, deseq.names.get('padj', 1), deseq.names.get('padjft', 1)
        )
        counts.compute(
          0, deseq.names.get('padj', 2), deseq.names.get('padjft', 2)
        )

        counts_norm = Ribopip::ExpressionAnalysis::Deseq::AddCountsToResults.new(
          deseq.names.get('merged'), '_norm', deseq.names.get('log', 1)
        )
        counts_norm.compute(
          1, deseq.names.get('allft', 1), deseq.names.get('allftnorm', 1)
        )
        counts_norm.compute(
          0, deseq.names.get('allft', 2), deseq.names.get('allftnorm', 2)
        )
        counts_norm.compute(
          1, deseq.names.get('padjft', 1), deseq.names.get('padjftnorm', 1)
        )
        counts_norm.compute(
          0, deseq.names.get('padjft', 2), deseq.names.get('padjftnorm', 2)
        )
        counts_norm.write(deseq.names.get('mergednorm', 1))

        # add gene names to gene ids
        id2name.translate(
          1,
          "#{deseq.names.get('padjftnorm', 1)}",
          "#{deseq.names.get('padjnames', 1)}"
        )
        id2name.translate(
          0,
          "#{deseq.names.get('padjftnorm', 2)}",
          "#{deseq.names.get('padjnames', 2)}"
        )
      end

      ##
      # POSTDESEQ
      # Parse significant genes from each experiment of each condition
      fp_experiments.length.times do | i |
        break if deseq_results[fp_experiments.length + i].nil?
        fp_base = deseq_results[i].names.base
        mrna_base = deseq_results[fp_experiments.length + i].names.base
        basenames = [
          Ribopip::Filenames::Postdeseq.new(fp_base),
          Ribopip::Filenames::Postdeseq.new(mrna_base)
        ]

        # copy significant results
        basenames.each do |names|
          run_cmd("cp #{names.get('padjftnorm', 1)} #{names.get('padjcp', 1)}")
          run_cmd("cp #{names.get('padjftnorm', 2)} #{names.get('padjcp', 2)}")
        end

        # add non-significant genes to either FP and total, if they were
        # significant in the other experiment
        important = Ribopip::ExpressionAnalysis::Deseq::ImportantGenes.new(
          basenames
        )
        important.compute

        # add gene names to gene ids
        id2name = GeneID2Name.new(options[:annotation])
        basenames.each do |names|
          id2name.translate(
            1, "#{names.get('padjcp', 1)}", "#{names.get('padjnames', 1)}"
          )
          Ribopip::ExpressionAnalysis::Deseq::Deseq2Excel.convert(
            "#{names.get('padjnames', 1)}", "#{names.get('padjxml', 1)}"
          )
          id2name.translate(
            0, "#{names.get('padjcp', 2)}", "#{names.get('padjnames', 2)}"
          )
          Ribopip::ExpressionAnalysis::Deseq::Deseq2Excel.convert(
            "#{names.get('padjnames', 2)}", "#{names.get('padjxml', 2)}"
          )
        end
      end

      ##
      # Detecting TE changes
      fp_experiments.length.times do |i|
        break if deseq_results[fp_experiments.length + i].nil?
        fp = deseq_results[i]
        mrna = deseq_results[fp_experiments.length + i]
        names_te = Ribopip::Filenames::TE.new(fp.names.base)
        Ribopip::ExpressionAnalysis::Deseq::TE.compute(
          fp.names.get('mergednorm'),
          mrna.names.get('mergednorm'),
          names_te.get('counts')
        )
        deseq_path = File.expand_path(
            '../../../scripts/run_deseq.R', __FILE__
        )

        cmd = "#{deseq_path} 1 #{names_te.get('counts')} 128" \
              " #{fp.conditions[0]} #{fp.conditions[1]} #{fp.reps_num.join(' ')} "\
              " #{names_te.get('all')} #{names_te.get('padj')} " \
              " > #{names_te.get('log')}"
        p "RUN TE: #{cmd}"
        `#{cmd}`

        counts = Ribopip::ExpressionAnalysis::Deseq::AddCountsToResults.new(
          names_te.get('counts'), '_TE'
        )
        counts.compute(1, names_te.get('padj'), names_te.get('padjft'))

        countsfp = Ribopip::ExpressionAnalysis::Deseq::AddCountsToResults.new(
          fp.names.get('merged'), '_fp'
        )
        countsfp.compute(1, names_te.get('padjft'), names_te.get('padjftfp'))

        countsmrna = Ribopip::ExpressionAnalysis::Deseq::AddCountsToResults.new(
          mrna.names.get('merged'), '_mrna'
        )
        countsmrna.compute(
          1, names_te.get('padjftfp'), names_te.get('padjftmrna')
        )

        id2name = GeneID2Name.new(options[:annotation])
        # add gene names to gene ids
        id2name.translate(
          1,
          names_te.get('padjftmrna'),
          names_te.get('padjnames')
        )

        Ribopip::ExpressionAnalysis::Deseq::Deseq2Excel.convert(
          names_te.get('padjnames'), names_te.get('xml')
        )
      end

      ##
      # RIBODIFF and DEseq with multiple factors
      # Run RiboDiff for each condition against each other condition
      fp_experiments.length.times do | i |
        (i+1).upto(fp_experiments.length-1) do |j|
          fp_cond1 = fp_experiments[i]
          fp_cond2 = fp_experiments[j]
          total_cond1 = mrna_experiments[i]
          total_cond2 = mrna_experiments[j]

          # RIBODIFF
          names_rdiff = Ribopip::Filenames::Ribodiff.new(
            fp_cond1[0].gsub(/fp_/, ''),
            fp_cond2[0].gsub(/fp_/, ''),
            options[:outdir]
          )
          rdiff = Ribopip::ExpressionAnalysis::Ribodiff.new(
            names_rdiff.dup,
            [fp_cond1[1..-1], fp_cond2[1..-1], total_cond1[1..-1], total_cond2[1..-1]],
            [fp_cond1[0], fp_cond2[0], total_cond1[0], total_cond2[0]],
            options[:outdir]
          )
          rdiff.compute_exp_outline
          rdiff.merge_featurecounts
          if options[:ribodiff]
            rdiff.compute
            counts = Ribopip::ExpressionAnalysis::Deseq::AddCountsToResults.new(
              names_rdiff.get('merged')
            )
            counts.compute(
              0, names_rdiff.get('padj'), names_rdiff.get('padjft')
            )
            id2name = GeneID2Name.new(options[:annotation])
            # add gene names to gene ids
            id2name.translate(
              0,
              names_rdiff.get('padjft'),
              names_rdiff.get('padjnames')
            )

            Ribopip::ExpressionAnalysis::Deseq::Deseq2Excel.convert(
              names_rdiff.get('padjnames'), names_rdiff.get('xml')
            )
          end

          # DESEQ WITH MULTI-FACTORS
          names_de = Ribopip::Filenames::Deseq.new(
            fp_cond1[0].gsub(/fp_/, ''),
            fp_cond2[0].gsub(/fp_/, ''),
            options[:outdir]
          )
          deseq_path = File.expand_path(
              '../../../scripts/run_deseq_rp.R', __FILE__
          )

          cmd = "#{deseq_path} #{names_rdiff.get('merged')}" \
                " #{names_rdiff.get('outline')} 128"\
                " #{names_de.get('all', '_mult')} "\
                " #{names_de.get('padj', '_mult')} " \
                " > #{names_de.get('log', '_mult')}"
          p "RUN DEseq RP: #{cmd}"
          `#{cmd}`

          counts = Ribopip::ExpressionAnalysis::Deseq::AddCountsToResults.new(
            names_rdiff.get('merged', '_mult')
          )
          counts.compute(
            0, names_de.get('padj', '_mult'), names_de.get('padjft', '_mult')
          )

          id2name = GeneID2Name.new(options[:annotation])
          # add gene names to gene ids
          id2name.translate(
            0,
            names_de.get('padjft', '_mult'),
            names_de.get('padjnames', '_mult')
          )

          Ribopip::ExpressionAnalysis::Deseq::Deseq2Excel.convert(
            names_de.get('padjnames', '_mult'), names_de.get('xml', '_mult')
          )
        end
      end

    end
  end
end
