module Ribopip
  # FeatureCounts
  class FeatureCounts
    extend Ribopip

    # Detects key collisions in hash and push elements into an array to avoid
    # data loss
    #
    # hash - hash
    # key  - key
    # data - data
    def self.push_into_hash(hash, key, data)
      if hash[key.to_sym].nil?
        hash[key.to_sym] = [data]
      else # existing key
        hash[key.to_sym].push(data)
      end
    end

    # Parses GTF annotation and writes a tab-seperated file, which assigns
    # biotypes to gene_ids. This file will be used to process FeatureCounts
    # output. Writes "#{annotation}.types".
    #
    # annotation - input gtf file
    def self.write_index(annotation)
      features = {}
      names = {}
      # get 9th column of all exons
      exons_annotation = `cat #{annotation} | grep -P '\texon\t' | cut -f9`

      # build hash
      exons_annotation.each_line do |line|
        gene_id = line[/gene_id\W+\w+/].delete('"').split[1]
        name = line[/gene_name\W+\w+/].delete('"').split[1]
        name = name.nil? ? 'n/a' : name
        type = line[/gene_biotype\W+\w+/].delete('"').split[1]

        # new gene_id
        push_into_hash(features, gene_id, type)
        push_into_hash(names, gene_id, name)
      end

      # write index
      file_types = File.open("#{annotation}.types", 'w')
      features.each do |gene_id, ft|
        file_types.printf("%s\t%s\t%s\n", gene_id, ft.uniq.sort * '; ',
                          names[:"#{gene_id}"].uniq.sort * '; ')
      end
      file_types.close
    end

    # Parses a tab-seperated file into a hash, which assigns biotypes to
    # gene_ids. This file will be used to process FeatureCounts output.
    #
    # filename - input file
    #
    # Returns hash.
    def self.read_index(filename)
      features = {}
      File.readlines(filename).each do |line|
        key = line.split("\t")[0]
        values = line.split("\t")[1].strip.split('; ')
        features[:"#{key}"] = values
      end

      features
    end

    # Computes feature distribution from feature counts output
    #
    # infile     - input file
    # outdistri  - output file containing ft distribution
    # outsums    - output file containing ft distribution sums
    # annotation - genomic annotation
    def self.compute_distri(infile, outdistri, outsums, annotation)
      fail "#{infile} does not exist." unless File.exist?(infile)

      # build feature index for annotation, if non-existent
      write_index(annotation) unless File.exist?("#{annotation}.types")

      # read feature index into hash
      features = read_index("#{annotation}.types")
      cumulative_sums = {}

      # write feature distribution
      file_distri = File.open(outdistri, 'w')
      File.readlines(infile).drop(1).each do |line|
        gene_id = line.split[0]
        type = features[:"#{gene_id}"].nil? ? 'nil' : features[:"#{gene_id}"] * '; '
        num = line.split[6].to_i
        next if num < 0
        file_distri.printf("%s\t%s\t%d\n", gene_id, type, num)
        if cumulative_sums[:"#{type}"].nil?
          cumulative_sums[:"#{type}"] = num
        else
          cumulative_sums[:"#{type}"] += num
        end
      end
      file_distri.close

      # write cumulative sums in descending order
      file_sums = File.open(outsums, 'w')
      cumulative_sums.sort_by(&:last).reverse_each do |type, num|
        file_sums.printf("%s\t%d\n", type, num)
      end
      file_sums.close
    end

    # Wrapper function. Runs feature counts and computes field distribution
    #
    # infile     - input file
    # outfile    - output file
    # annotation - genomic annotation (gtf)
    # target     - gtf feature type (CDS, exon, gene, Selenocysteine,
    #              start_codon, stop_codon, transcript, UTR)
    def self.compute(infile, outfile, annotation, target = 'exon')
      fail "#{infile} does not exist." unless File.exist?(infile)
      run_cmd(
        "featureCounts -t #{target} -a #{annotation} -o #{outfile} #{infile}"
      )
    end
  end
end
