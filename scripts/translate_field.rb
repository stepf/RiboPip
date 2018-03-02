#!/usr/bin/env ruby

usespace_d=ARGV.delete("-sd")
usespace_i=ARGV.delete("-si")
if ARGV.delete("-s")
  usespace_d=true
  usespace_i=true
end

def parse_argv
  args = {}
  if ARGV.size != 3
    # %-PURPOSE-%
    purpose = "Translate a column in a TSV file using a dictionary"
    STDERR.puts purpose
    STDERR.puts "Usage: #$0 [-s[d|i][-n] <input_tsv_file> <acc_dictionary> "+
                           "<1-based-column_to_convert>"
    STDERR.puts "-sd dictionary uses space as delimiter instead of tabs"
    STDERR.puts "-si input uses space as delimiter instead of tabs"
    STDERR.puts "-s equivalent to -sd -si"
    exit 1
  end
  args[:input_filename]=ARGV.shift
  args[:dictionary_filename]=ARGV.shift
  args[:column]=Integer(ARGV.shift)-1
  return args
end

def prepare_dictionary(dict_filename, usespace)
  d=File.open(dict_filename)
  dictionary={}
  d.each do |line|
    if usespace
      elems=line.split(" ")
      k=elems.shift
      elems[-1].chomp!
      dictionary[k]=elems
    else
      elems=line.split("\t")
      dictionary[elems[0].to_sym]=elems[1].chomp
    end
  end
  d.close
  return dictionary
end

def apply_dictionary(input_filename, column, dictionary, usespace)
  f=File.open(input_filename)
  f.each do |line|
    elems=line.split(usespace ? " " : "\t")
    translated=dictionary[elems[column].chomp.to_sym]
    elems[column]=translated if !translated.nil?
    puts elems.join(usespace ? " " : "\t")
  end
end

args = parse_argv
apply_dictionary(args[:input_filename],
                 args[:column],
                 prepare_dictionary(args[:dictionary_filename],
                                    usespace_d),
                 usespace_i)
