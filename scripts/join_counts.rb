#!/usr/bin/env ruby

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

# write tsv
File.open(tsv_out, 'w')
tsv_out.puts header.join("\t")
data.each do |elems|
  tsv_out.puts elems.join("\t")
end
tsv_out.close

# write tsv
File.open(xml_out, 'w')
xml_out.puts('<?xml version="1.0"?>')
xml_out.puts('<Workbook xmlns="urn:schemas-microsoft-com:office:' \
              'spreadsheet"')
xml_out.puts('  xmlns:o="urn:schemas-microsoft-com:office:office"')
xml_out.puts('  xmlns:x="urn:schemas-microsoft-com:office:excel"')
xml_out.puts('  xmlns:ss="urn:schemas-microsoft-com:office:spreadsheet"')
xml_out.puts('  xmlns:html="http://www.w3.org/TR/REC-html40">')
xml_out.puts('  <Worksheet ss:Name="Sheet1">')
xml_out.puts('    <Table>')

xml_out.puts('      <Row>')
header.each do |e|
  t = is_float?(e) ? 'Number' : 'String'
  xml_out.puts("        <Cell><Data ss:Type=\"#{t}\">#{e}</Data></Cell>")
end
xml_out.puts('      </Row>')

data.each do |row|
  xml_out.puts('      <Row>')
  row.each do |e|
    t = is_float?(e) ? 'Number' : 'String'
    xml_out.puts("        <Cell><Data ss:Type=\"#{t}\">#{e}</Data></Cell>")
  end
  xml_out.puts('      </Row>')
end

xml_out.puts('    </Table>')
xml_out.puts('  </Worksheet>')
xml_out.puts('</Workbook>')
xml_out.close
