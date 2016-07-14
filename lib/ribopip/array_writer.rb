module Ribopip
  # Writes arrays into the following formats:
  #   - TSV
  #   - ExcelXML
  module ArrayWriter
    # abstract writer class
    class Writer
      # Initializes duplicate of input array to prevent unwanted side effects
      #
      # array         - input array
      # is_countarray - boolean
      def initialize(array, is_countarray = false)
        @array = array.dup
        @array.insert(0, ['Name', '#Reads', 'Perc']) if is_countarray
      end

      # write array to output file
      #
      # outfile - output filename
      def write(outfile)
        file_out = File.open(outfile, 'w')

        write_header(file_out)
        write_content(file_out)
        write_footer(file_out)

        file_out.close
      end

      # abstract method
      def write_header(_)
      end

      # Writes array to file
      #
      # file_out - file object
      def write_content(file_out)
        file_out.puts(@array)
      end

      # abstract method
      def write_footer(_)
      end

      # checks if string contains a float
      #
      # str - input string
      #
      # Returns boolean
      def float?(str)
        true if Float(str) rescue false
      end
    end

    # concrete class for tab-seperated files
    class TSVWriter < Writer
      def write_content(file_out)
        @array.each do |row|
          row.each do |elem|
            elem = elem.to_f.round(2) if float?(elem)
            file_out.printf("%s\t", elem)
          end
          file_out.printf("\n")
        end
      end
    end

    # concrete class for ExcelXML, for specification see:
    # https://en.wikipedia.org/wiki/Microsoft_Office_XML_formats
    class XMLWriter < Writer
      # begin workbook and a single sheet
      #
      # file_out - output file (class: File)
      #
      # Returns nothing
      def write_header(file_out)
        file_out.puts('<?xml version="1.0"?>')
        file_out.puts('<Workbook xmlns="urn:schemas-microsoft-com:office:' \
                      'spreadsheet"')
        file_out.puts('  xmlns:o="urn:schemas-microsoft-com:office:office"')
        file_out.puts('  xmlns:x="urn:schemas-microsoft-com:office:excel"')
        file_out.puts(
          '  xmlns:ss="urn:schemas-microsoft-com:office:spreadsheet"'
        )
        file_out.puts('  xmlns:html="http://www.w3.org/TR/REC-html40">')
        file_out.puts('  <Worksheet ss:Name="Sheet1">')
        file_out.puts('    <Table>')
      end

      # Iterates over input_array and writes rows / cells
      #
      # file_out - output file (class: File)
      def write_content(file_out)
        @array.each do |row|
          file_out.puts('      <Row>')
          row.each do |e|
            t = float?(e) ? 'Number' : 'String'
            file_out.puts(
              "        <Cell><Data ss:Type=\"#{t}\">#{e}</Data></Cell>"
            )
          end
          file_out.puts('      </Row>')
        end
      end

      # Closes every open tag
      #
      # file_out - output file (class: File)
      def write_footer(file_out)
        file_out.puts('    </Table>')
        file_out.puts('  </Worksheet>')
        file_out.puts('</Workbook>')
      end
    end
  end
end
