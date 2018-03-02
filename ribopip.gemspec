## This is the rakegem gemspec template. Make sure you read and understand
## all of the comments. Some sections require modification, and others can
## be deleted if you don't need them. Once you understand the contents of
## this file, feel free to delete any comments that begin with two hash marks.
## You can find comprehensive Gem::Specification documentation, at
## http://docs.rubygems.org/read/chapter/20
Gem::Specification.new do |s|
  s.specification_version = 2 if s.respond_to? :specification_version=
  s.required_rubygems_version = Gem::Requirement.new('>= 0') if s.respond_to? :required_rubygems_version=
  s.rubygems_version = '1.3.5'

  ## Leave these as is they will be modified for you by the rake gemspec task.
  ## If your rubyforge_project name is different, then edit it and comment out
  ## the sub! line in the Rakefile
  s.name              = 'ribopip'
  s.version           = '0.1.0'
  s.date              = '2017-12-15'
  s.rubyforge_project = 'ribopip'

  ## Make sure your summary is short. The description may be as long
  ## as you like.
  s.summary     = 'An alignment and analysis pipeline for Ribosome Profiling' \
                  'Data.'
  s.description = 'Long description. Maybe copied from the README.'

  ## List the primary authors. If there are a bunch of authors, it's probably
  ## better to set the email to an email list or something. If you don't have
  ## a custom homepage, consider using your GitHub URL or the like.
  s.authors  = ['Stefan Dang']
  s.email    = 'stefan@man-dangt.de'
  s.homepage = 'http://github.com/stepf'

  ## This gets added to the $LOAD_PATH so that 'lib/NAME.rb' can be required as
  ## require 'NAME.rb' or'/lib/NAME/file.rb' can be as require 'NAME/file.rb'
  s.require_paths = %w(lib)

  ## This sections is only necessary if you have C extensions.
  # s.require_paths << 'ext'
  # s.extensions = %w(ext/fastq_bucketize/extconf.rb)

  ## If your gem includes any executables, list them here.
  # s.executables = ["name"]

  ## Specify any RDoc options here. You'll want to add your README and
  ## LICENSE files to the extra_rdoc_files list.
  s.rdoc_options = ['--charset=UTF-8']
  #s.extra_rdoc_files = %w[README LICENSE]

  ## List your runtime dependencies here. Runtime dependencies are those
  ## that are needed for an end user to actually USE your code.
  s.add_dependency('thor', ['>= 0.19.1'])

  ## List your development dependencies here. Development dependencies are
  ## those that are only needed during development
  # s.add_development_dependency('DEVDEPNAME', [">= 1.1.0", "< 2.0.0"])
  s.add_development_dependency('rdoc', ['>= 3.2.3'])
  s.add_development_dependency('rspec', ['>= 3.5.0'])
  s.add_development_dependency('simplecov', ['>= 0.12.0'])

  ## Leave this section as-is. It will be automatically generated from the
  ## contents of your Git repository via the gemspec task. DO NOT REMOVE
  ## THE MANIFEST COMMENTS, they are used as delimiters by the task.
  # = MANIFEST =
  s.files = %w[
    Dockerfile
    Gemfile
    LICENSE
    Rakefile
    Readme.md
    bin/ribopip
    ext/fastq-bucketize-0.1/Makefile
    ext/fastq-bucketize-0.1/fastq-bucketize.c
    lib/ribopip.rb
    lib/ribopip/array_writer.rb
    lib/ribopip/binary_checker.rb
    lib/ribopip/cli.rb
    lib/ribopip/counts.rb
    lib/ribopip/deseq.rb
    lib/ribopip/expression_analysis.rb
    lib/ribopip/feature_counts.rb
    lib/ribopip/filenames.rb
    lib/ribopip/gene_id_2_name.rb
    lib/ribopip/metrics.rb
    lib/ribopip/pipeline.rb
    ribopip.gemspec
    scripts/bootstrap
    scripts/join_counts.rb
    scripts/run_deseq.R
    scripts/translate_field.rb
    spec/array_writer_spec.rb
    spec/binary_checker_spec.rb
    spec/counts_spec.rb
    spec/expression_analysis_spec.rb
    spec/filenames_spec.rb
    spec/gene_id_2_name_spec.rb
    spec/metrics_spec.rb
    spec/pipeline_spec.rb
    spec/ribseq_spec.rb
    spec/spec_helper.rb
    spec/testdata/SRR315601_625.fastq
    spec/testdata/array_writer/SRR315601_625.counts.tsv
    spec/testdata/array_writer/SRR315601_625.counts.txt
    spec/testdata/array_writer/SRR315601_625.counts.xml
    spec/testdata/array_writer/expected_read_counts.tsv
    spec/testdata/array_writer/expected_read_counts.txt
    spec/testdata/array_writer/expected_read_counts.xml
    spec/testdata/counts/align_summary.txt
    spec/testdata/counts/uni_mapped_0err.bam
    spec/testdata/expression_analysis/d1_vs_d2.counts
    spec/testdata/expression_analysis/d1_vs_d2.deseq.counts.norm.tsv
    spec/testdata/expression_analysis/d1_vs_d2.deseq.counts.tsv
    spec/testdata/expression_analysis/d1_vs_d2.deseq1.all.counts.csv
    spec/testdata/expression_analysis/d1_vs_d2.deseq1.all.counts.norm.csv
    spec/testdata/expression_analysis/d1_vs_d2.deseq1.all.csv
    spec/testdata/expression_analysis/d1_vs_d2.deseq1.all.ft.csv
    spec/testdata/expression_analysis/d1_vs_d2.deseq1.all.ft.norm.csv
    spec/testdata/expression_analysis/d1_vs_d2.deseq1.log
    spec/testdata/expression_analysis/d1_vs_d2.deseq1.padj.csv
    spec/testdata/expression_analysis/d1_vs_d2.deseq1.padj.ft.csv
    spec/testdata/expression_analysis/d1_vs_d2.deseq2.log
    spec/testdata/expression_analysis/d1_vs_d2.rdiff.outline.csv
    spec/testdata/expression_analysis/d1_vs_d2.rdiff.results.padj.tsv
    spec/testdata/expression_analysis/d1_vs_d2.rdiff.results.tsv
    spec/testdata/expression_analysis/dummy1_rep1.txt
    spec/testdata/expression_analysis/dummy1_rep2.txt
    spec/testdata/expression_analysis/dummy1_rep3.txt
    spec/testdata/expression_analysis/dummy2_rep1.txt
    spec/testdata/expression_analysis/dummy2_rep2.txt
    spec/testdata/expression_analysis/expected.all.ft.csv
    spec/testdata/expression_analysis/expected.all.ft.norm.csv
    spec/testdata/expression_analysis/expected.counts.norm.tsv
    spec/testdata/expression_analysis/expected.rdiff.outline.csv
    spec/testdata/expression_analysis/expected.rdiff.results.padj.tsv
    spec/testdata/expression_analysis/expected_d1_vs_d2.counts
    spec/testdata/expression_analysis/expected_names.xml
    spec/testdata/expression_analysis/mrna.deseq1.all.ft.csv
    spec/testdata/expression_analysis/mrna.deseq1.padj.ft.csv
    spec/testdata/expression_analysis/names.csv
    spec/testdata/expression_analysis/names.xml
    spec/testdata/gene_id_2_name/d1_vs_d2.deseq1.all.ft.csv
    spec/testdata/gene_id_2_name/d1_vs_d2.deseq1.names.csv
    spec/testdata/gene_id_2_name/expected.names.csv
    spec/testdata/ref/expected.name_index
    spec/testdata/ref/mini.GRCm38.80.gtf
    spec/testdata/ref/mini.GRCm38.80.gtf.name_index
    spec/testdata/statistics/SRR315601_625.fastq.gt_seqstat
    spec/testdata/statistics/expected.gt_seqstat
  ]
  # = MANIFEST =

  ## Test files will be grabbed from the file list. Make sure the path glob
  ## matches what you actually use.
  s.test_files = s.files.select { |path| path =~ /^test\/test_.*\.rb/ }
end
