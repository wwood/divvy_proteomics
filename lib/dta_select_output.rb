




module Bio::DTASelect

  class OutputFile
    def self.log
      SelectedProtein.new.log
    end

    class SelectedProtein
      include Bio::DivvyProteomics::Logging
      include Bio::DivvyProteomics::DivvyableProtein

      attr_accessor :identifier

      attr_accessor :sequence_count, :spectrum_count, :sequence_coverage, :length, :molwt, :pi, :validation_status, :descriptive_name

      attr_accessor :peptides

      def initialize
        @peptides = []
      end



      def log
        Bio::Log::LoggerPlus[LOG_NAME]
      end
    end

    class Peptide
      include Bio::DivvyProteomics::Logging

      attr_accessor :identifier

      # Hash of column names to values. These are different for different DTAselect output files, it seems.
      attr_accessor :dtaselect_attributes

      # Array of proteins that have this peptide associated
      attr_accessor :parent_proteins

      def initialize
        @parent_proteins = []
      end

      def inspect
        "Peptide: #{@parent_proteins.length} @parent_proteins: [#{@parent_proteins.collect{|pro| pro.identifier}.join(', ')} @identifier: #{identifier}, @attributes: #{dtaselect_attributes.inspect}]"
      end

      def redundancy
        @dtaselect_attributes['Redundancy'].to_i
      end

      def reported_unique?
        dtaselect_attributes.length == 1
      end
    end

    class Result
      include Bio::DivvyProteomics::Logging

      # hash of protein identifier to Protein object
      attr_accessor :protein_name_to_object

      # hash of peptide identifier to Peptide object
      attr_accessor :peptide_name_to_object
    end

    def self.parse(io)
      result = Result.new

      # Hashes of identifiers to objects
      result.protein_name_to_object = {}
      result.peptide_name_to_object = {}

      # Read in the tab separated file
      reading_header = true
      current_proteins = []
      last_line_was_protein_name = false
      peptide_attribute_names = nil

      # Parse each line of the DTAselect file
      io.each_line do |line|
        splits = line.chomp.split("\t")
        log.debug "Parsing line `#{line.chomp}'" if log.debug?

        if reading_header
          log.debug "reading header"
          if splits[0] == 'Unique'
            reading_header = false

            # Current line describes the peptide attributes
            peptide_attribute_names = splits

            # This field has special importance, so be picky
            raise "Badly parsed file at this line: #{line.inspect}, expected 2nd field to be 'FileName', found #{splits[1]}" unless splits[1] == 'FileName'
          end
          next
        end

        # OK, now we are reading the actual table, not the header
        if splits[0] != '' and splits[11].nil?
          ident = splits[0]

          if !last_line_was_protein_name
            # Sometimes several proteins are given all in the one header line
            # start a new protein
            log.debug "New protein now being parsed" if log.debug?
            current_proteins = []
          end

          current_protein = SelectedProtein.new
          last_line_was_protein_name = true
          current_proteins.push current_protein

          current_protein.identifier = ident

          i = 1
          current_protein.sequence_count = splits[i].to_i; i+=1
          current_protein.spectrum_count = splits[i].to_i; i+=1
          current_protein.sequence_coverage = splits[i].to_f; i+=1
          current_protein.length = splits[i].to_i; i+=1
          current_protein.molwt = splits[i].to_f; i+=1
          current_protein.pi = splits[i].to_f; i+=1
          current_protein.validation_status = splits[i].to_f; i+=1
          current_protein.descriptive_name = splits[i]

          if result.protein_name_to_object[ident]
            raise "Unexpectedly found the same protein identifier twice: #{ident}, from line #{line.chomp}"
          end
          result.protein_name_to_object[ident] = current_protein



        elsif splits[1] == 'Proteins'
          # Done processing, except for the bits down the bottom which aren't parsed (yet, at least)
          break



        else
          log.debug "New spectra now being parsed" if log.debug?
          last_line_was_protein_name = false

          # Record a spectra
          ident = splits[1]
          raise "Unexpected hits name `#{ident}', from line `#{line.chomp}'" unless ident.length > 10

          pep = result.peptide_name_to_object[ident]
          if pep.nil?
            pep = Peptide.new
            pep.identifier = ident

            peptide_attribute_names.each_with_index do |attribute_name,i|
              pep.dtaselect_attributes ||= {}
              pep.dtaselect_attributes[attribute_name] = splits[i]
            end

            result.peptide_name_to_object[ident] = pep
          end

          current_proteins.each do |current_protein|
            pep.parent_proteins.push current_protein
            current_protein.peptides.push pep
          end
          log.debug "Parsed this peptide #{pep.inspect}" if log.debug?
        end
      end

      log.debug "Proteins parsed: #{result.protein_name_to_object.inspect}" if log.debug?
      return result
    end
  end
end
