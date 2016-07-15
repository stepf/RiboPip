module Ribopip
  # Parses number of reads from various files and stores them in an count array
  module Counts
    # keep track of #reads using the Counts::CountArray and
    # Counts::Parser classes
    class Bookkeeper
      # array - array; read counts
      attr_reader :array

      # Initiliazes empty Counts::Array
      def initialize
        @array = Ribopip::Counts::Array.new
      end

      # parses #reads from file, pushes into count array
      #
      # infile   - input filename
      # name     - name of #reads
      # regex    - optional: regex pattern
      #
      # Returns nothing
      def parse_nreads(infile, name, regex = nil)
        fail "#{infile} does not exist." unless File.exist?(infile)
        filetype = regex.nil? ? File.extname(infile)[1..-1] : 'txt'

        # call lambda according to filetype to get #reads
        lambda = { sam: FLAGSTAT, bam: FLAGSTAT, fastq: LINECOUNT, txt: REGEX }
        num = lambda[filetype.to_sym].call(infile, regex)

        @array.push(name, num)
      end
    end

    # defining the CountArray:
    # [[name1, count1, perc1], [name2, count2, perc2], ...]
    class Array
      # the count array
      attr_reader :array

      # initializes (empty) array
      def initialize
        @array = []
      end

      # overwriting Array.each method
      def each
        @array.each do |tuple|
          fail 'CountArray seems to be corrupt' if tuple.length != 3
          yield tuple
        end
      end

      # overwriting Array.insert method
      def insert(*args)
        @array.insert(*args)
      end

      # overwriting Array.dup method
      def dup
        @array.dup
      end

      # pushes #reads into count array; automatically computes percentage w.r.t.
      # first entry of the array
      #
      # name  - count name
      # count - #reads
      def push(name, count)
        perc = @array.empty? ? 100 : count.to_f / @array[0][1] * 100
        @array.push([name, count, perc])
      end
    end

    # parsing #reads from files
    class Parser
      # get #reads using samtools flagstat for sam / bam files
      #
      # infile - input file
      #
      # Returns integer
      def self.flagstat(infile)
        `samtools flagstat #{infile}`[/^.*mapped/].split.first.to_i
      end

      # get #reads using linecount for fastq
      #
      # infile - input file
      #
      # Returns integer
      def self.linecount(infile)
        `wc -l #{infile}`.split.first.to_i / 4
      end

      # get #reads using regex pattern
      #
      # infile - input file
      # regex  - regex pattern
      #
      # Returns integer
      def self.regex(infile, regex)
        number = File.read(infile).delete(',')[/#{regex}.*$/]
        number.nil? ? 0 : number[/\d+/].to_i
      end
    end

    # get #reads using samtools flagstat for sam / bam files
    #
    # infile - input file
    #
    # Returns integer
    FLAGSTAT = lambda do |infile, _|
      `samtools flagstat #{infile}`[/^.*mapped/].split.first.to_i
    end

    # get #reads using linecount for fastq
    #
    # infile - input file
    #
    # Returns integer
    LINECOUNT = lambda do |infile, _|
      `wc -l #{infile}`.split.first.to_i / 4
    end

    # get #reads using regex pattern
    #
    # infile - input file
    # regex  - regex pattern
    #
    # Returns integer
    REGEX = lambda do |infile, regex|
      number = File.read(infile).delete(',')[/#{regex}.*$/]
      number.nil? ? 0 : number[/\d+/].to_i
    end
  end
end
