module Ribopip
  # Loads gene names from an GTF annotation and adds them to any desired
  # delimited text file.
  class GeneID2Name
    # id2name - hash with key-value pairs of gene IDs to gene names
    attr_reader :id2name

    # annot - path to gtf annotation
    def initialize(annot)
      @annot = annot
      @id2name = {}
      @index = "#{@annot}.name_index"
      if File.exist?(@index)
        read_index
      else
        build_index
      end
    end

    # Builds @id2name hash and writes @index file for reusing it
    def build_index
      # parse annotation
      File.readlines(@annot).each do |line|
        next if line[0..1] == '#!'
        gene_id = line[/gene_id\W+\w+/].delete('"').split[1]
        gene_name = line[/gene_name\W+\w+/].delete('"').split[1]
        next if gene_name.nil?
        @id2name[gene_id.to_sym] = gene_name
      end

      # write index
      index = File.open(@index, 'w')
      @id2name.each do |gene_id, gene_name|
        index.printf("%s\t%s\n", gene_id, gene_name)
      end
      index.close
    end

    # Reads pre-built index
    def read_index
      File.readlines(@index).each do |line|
        gene_id, gene_name = line.chomp.split("\t")
        @id2name[gene_id.to_sym] = gene_name
      end
    end

    # Translates gene_ids in tab-seperated files to gene_names and writes an
    # additional column into a new file.
    #
    # id_col   - column containing the ID, starting with 0
    # infile   - input file
    # outfile  - output file
    def translate(id_col, infile, outfile, delimiter = "\t")
      out = File.open(outfile, 'w')
      File.readlines(infile).each do |line|
        inelems = line.split(delimiter)
        gene_id = inelems[id_col]
        if gene_id.nil?
          out.puts
          next
        end
        if !gene_id.empty?
          gene_id.delete!('"')
          gene_name = @id2name[gene_id.to_sym].nil? ? ' ' : @id2name[gene_id.to_sym]
        else # not available
          gene_name = ''
        end
        inelems[id_col] = [gene_id, gene_name]
        out.puts inelems.flatten!.join(delimiter)
      end
      out.close
    end
  end
end
