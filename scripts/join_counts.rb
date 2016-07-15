#!/usr/bin/env ruby

require_relative '../lib/ribopip/array_writer.rb'

if (ARGV.size == 0) or (ARGV.size % 2 != 0)
  puts "Usage: #{$0} <file1> <label1> <file2> <label2> [<file3> <label3>...] <tsv out> <xml out>"
  exit 1
end
files = []
labels = []
until ARGV.empty?
  if ARGV.size == 2
    tsv_out = ARGV.shift
    xml_out = ARGV.shift
  end
  files << ARGV.shift
  labels << ARGV.shift
end
n_spaces = []
data = []
files.each do |filename|
  n_spaces_determined = false
  File.readlines(filename).each do |line|
    elems = line.chomp.split("\t")

    unless n_spaces_determined
      n_spaces << (elems.size - 2)
      n_spaces_determined = true
    end

    data << [elems[0]] if data.assoc(elems[0]).nil?

    elems[1..-1].each do |elem|
      data.assoc(elems[0]) << elem
    end
  end
end
header = ['']
i = 0
until n_spaces.empty?
  header << labels[i]
  n_spaces.shift.times do
    header << ''
  end
  i += 1
end

outarray = [header] + data
tsv = Ribopip::ArrayWriter::TSVWriter.new(outarray)
tsv.write(tsv_out)
xml = Ribopip::ArrayWriter::XMLWriter.new(outarray)
xml.write(xml_out)
