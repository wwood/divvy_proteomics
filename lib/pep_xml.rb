require 'rexml/document'

class Bio::PepXML
  include Bio::DivvyProteomics::Logging

  attr_accessor :protein_name_to_object, :peptide_name_to_object

  class Protein
    include Bio::DivvyProteomics::Logging
    include Bio::DivvyProteomics::DivvyableProtein

    # Array of peptide objects that have been assigned to this protein
    attr_accessor :peptides

    attr_accessor :identifier, :descriptive_name
  end

  # Named 'Peptide' but really mean Spectra. Just too hard to change
  class Peptide
    attr_accessor :parent_proteins

    # Name of the spectra
    attr_accessor :identifier

    def initialize
      @parent_proteins = []
    end

    #TODO: right now this just always returns 1. It should really be working out redundancy
    #properly by comparison of peptide sequences, but this isn't yet parsed this info
    def redundancy
      1
    end
  end

  def self.log
    Bio::PepXML.new.log
  end

  def self.parse(io)
    protein_name_to_object = {}
    peptide_name_to_object = {}

    #pep.elements.each('msms_pipeline_analysis/msms_run_summary/spectrum_query/search_result/search_hit'){|e|
    #  c+=1; p e.attributes['protein_descr'].strip;
    #  e.elements.each{|e|
    #    p e.name, e.attributes['protein_descr'].strip};break}
    xml = REXML::Document.new(io)

    parse_name_and_description = lambda do |e|
      name = e.attributes['protein'].strip
      description = e.attributes['protein_descr'].strip
      if name.nil? or name == ''
        name = e.attributes['protein_descr'].strip
      else
        description = name+' '+description
      end
      name.gsub!(/\t.*/,'')
      description.gsub!(/[\t\n]/,' ')

      [name, description]
    end

    #TODO: some better sanity checking here would be ideal.
    num_hits_parsed = 0
    xml.elements.each('msms_pipeline_analysis/msms_run_summary/spectrum_query/search_result/search_hit') do |hit|
      hit_number = hit.attributes['hit_rank']
      raise "Parsing error on #{hit}" if hit_number.nil?
      next if hit_number != "1"

      # Parse the primary hit
      name1, description1 = parse_name_and_description.call(hit)
      raise "No protein name found in this xml fragment: #{hit.to_s}" if name1.nil?
      spectrum_name = hit.parent.parent.attributes['spectrum'].strip
      raise "Parsing error (couldn't find spectrum name) with spectra #{hit.inspect}" if spectrum_name.nil?

      # It is possible to have multiple peptides both hit the spectra with hit_rank="1"
      # This happens when when e.g. leucine and isoleucine are possible.
      spectrum = peptide_name_to_object[spectrum_name]
      if spectrum.nil?
        spectrum = Peptide.new
        spectrum.identifier = spectrum_name
        peptide_name_to_object[spectrum_name] = spectrum
      end


      protein1 = protein_name_to_object[name1]
      if protein1.nil?
        protein1 = Protein.new
        protein1.identifier = name1
        protein1.descriptive_name = description1
        protein1.peptides = []
        protein_name_to_object[name1] = protein1
      end
      protein1.peptides.push spectrum
      spectrum.parent_proteins ||= []
      spectrum.parent_proteins.push protein1


      # Parse the alternate hits. Only look at children with protein_descr attributes - these are
      # these are the alternate proteins
      hit.each_element_with_attribute('protein_descr') do |e|
        name, description = parse_name_and_description.call(e)

        alternate = protein_name_to_object[name]
        if alternate.nil?
          alternate = Protein.new
          alternate.identifier = name
          alternate.descriptive_name = description
          alternate.peptides = []
          protein_name_to_object[name] = alternate
        end
        alternate.peptides.push spectrum
        spectrum.parent_proteins.push alternate
      end

      # Don't count the same protein multiple times - might happen when a spectru
      spectrum.parent_proteins.uniq!

      num_hits_parsed += 1
    end
    log.info "Parsed #{num_hits_parsed} search hits"

    pepxml = Bio::PepXML.new
    pepxml.protein_name_to_object = protein_name_to_object
    pepxml.peptide_name_to_object = peptide_name_to_object

    return pepxml
  end
end
