require 'rexml/document'

class Bio::PepXML
  include Bio::DivvyProteomics::Logging

  attr_accessor :protein_name_to_object, :peptide_name_to_object

  class Protein
  end

  class Peptide
    attr_accessor :proteins

    def initialize
      @proteins = []
    end
  end

  def self.parse(io)
    @protein_name_to_object = {}
    @peptide_name_to_object = {}

    #pep.elements.each('msms_pipeline_analysis/msms_run_summary/spectrum_query/search_result/search_hit'){|e|
    #  c+=1; p e.attributes['protein_descr'].strip;
    #  e.elements.each{|e|
    #    p e.name, e.attributes['protein_descr'].strip};break}
    xml = REXML::Document.new(io)

    #TODO: some better sanity checking here would be ideal.
    num_hits_parsed = 0
    xml.each('msms_pipeline_analysis/msms_run_summary/spectrum_query/search_result/search_hit') do |hit|
      hit_number = hit.attributes['hit_rank']
      raise "Parsing error on #{hit}" if hit_number.nil?
      next if hit_number != "1"

      name1 = hit.attributes['protein_descr']
      raise "No protein name found in this xml fragment: #{hit.to_s}" if name1.nil?
      spectra = Spectra.new
      protein1 = @protein_name_to_object[name1]
      if protein1.nil?
        protein1 = Protein.new
        protein1.descriptive_name = name1
        protein1.peptides = []
      end


      # Parse the alternate hits


      num_hits_parsed += 1
    end
    log.info "Parsed #{num_hits_parsed} search hits"
  end
end
