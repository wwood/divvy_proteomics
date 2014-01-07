require 'systemu'
require 'pp'
require 'open3'
require 'tempfile'

require 'spec_helper'



describe 'pepxml parsing' do
  let(:header){"ID\tUnique spectra\tNon-unique spectra\tEstimated total spectra\tNormalised spectral count\tDescription\tProteins sharing spectra\n"}
  it 'should parse decently' do
    pepxml = Bio::PepXML.parse(File.open(File.join(TEST_DATA_DIR, 'minimal.pep.xml')))

#        <spectrum_query spectrum="Tara_38sfc_FASP_8hr_OrbiVelosPro_Run1_030513_02.9921.9921.2" start_scan="9921" end_scan="9921" retention_time_sec="41.732728333333334" activation_method="CID" precursor_intensity="61015.84375" precursor_neutral_mass="1246.6412829403125" assumed_charge="2" index="1">
#        <search_hit hit_rank="1" peptide="IADQTIGTANSR" protein="" num_tot_proteins="3" num_matched_ions="0" calc_neutral_pep_mass="1246.6412829403125" massdiff="0" protein_descr="&gt;38SUR_2379_1524213_2&#x9;" protein_mw="43.185399974660044" calc_pI="5.63037109375">
#          <alternative_protein protein="" protein_descr="&gt;38SUR_6350_1528184_1&#x9;" protein_mw="24.663561404659987" />
#          <alternative_protein protein="" protein_descr="&gt;38SUR_80622_1602456_1&#x9;" protein_mw="30.364007294659981" />
#          <search_score name="XCorr" value="4.7916374206542969" />

    pepxml.kind_of?(Bio::PepXML).should == true

    pepxml.protein_name_to_object.keys.sort.should == [
      '>38SUR_2379_1524213_2',
      '>38SUR_6350_1528184_1',
      '>38SUR_80622_1602456_1',
    ].sort
    pepxml.peptide_name_to_object.keys.sort.should == [
      'Tara_38sfc_FASP_8hr_OrbiVelosPro_Run1_030513_02.9921.9921.2'
    ]
    pepxml.protein_name_to_object.values.each do |prot|
      prot.kind_of?(Bio::PepXML::Protein).should == true
    end
    pepxml.peptide_name_to_object.values.each do |prot|
      prot.kind_of?(Bio::PepXML::Peptide).should == true
    end

    prot1 = pepxml.protein_name_to_object['>38SUR_2379_1524213_2']
    prot1.identifier.should == '>38SUR_2379_1524213_2'
    prot1.descriptive_name.should == '>38SUR_2379_1524213_2'
  end

  it 'should respond to divvy proteomics module things' do
    pepxml = Bio::PepXML.parse(File.open(File.join(TEST_DATA_DIR, 'minimal.pep.xml')))

#        <spectrum_query spectrum="Tara_38sfc_FASP_8hr_OrbiVelosPro_Run1_030513_02.9921.9921.2" start_scan="9921" end_scan="9921" retention_time_sec="41.732728333333334" activation_method="CID" precursor_intensity="61015.84375" precursor_neutral_mass="1246.6412829403125" assumed_charge="2" index="1">
#        <search_hit hit_rank="1" peptide="IADQTIGTANSR" protein="" num_tot_proteins="3" num_matched_ions="0" calc_neutral_pep_mass="1246.6412829403125" massdiff="0" protein_descr="&gt;38SUR_2379_1524213_2&#x9;" protein_mw="43.185399974660044" calc_pI="5.63037109375">
#          <alternative_protein protein="" protein_descr="&gt;38SUR_6350_1528184_1&#x9;" protein_mw="24.663561404659987" />
#          <alternative_protein protein="" protein_descr="&gt;38SUR_80622_1602456_1&#x9;" protein_mw="30.364007294659981" />
#          <search_score name="XCorr" value="4.7916374206542969" />

    pepxml.kind_of?(Bio::PepXML).should == true

    prot1 = pepxml.protein_name_to_object['>38SUR_2379_1524213_2']
    prot1.peptides.length.should == 1
    prot1.unique_spectra.should == 0
    prot1.non_unique_spectra.should == 1
    prot1.estimated_spectral_count.should == 0.0


    prot1 = pepxml.protein_name_to_object['>38SUR_6350_1528184_1']
    prot1.peptides.length.should == 1
    prot1.unique_spectra.should == 0
    prot1.non_unique_spectra.should == 1
    prot1.estimated_spectral_count.should == 0.0
  end

  it 'should respond to divvy proteomics module things with 1 unique hit' do
    pepxml = Bio::PepXML.parse(File.open(File.join(TEST_DATA_DIR, 'minimal2.pep.xml')))
    pepxml.kind_of?(Bio::PepXML).should == true

    prot1 = pepxml.protein_name_to_object['>38SUR_2379_1524213_2']
    prot1.peptides.length.should == 1
    prot1.unique_spectra.should == 1
    prot1.non_unique_spectra.should == 0
    prot1.estimated_spectral_count.should == 1.0
  end

  it 'should respond to divvy proteomics module things with 2 hits, where 1 is unique' do
    pepxml = Bio::PepXML.parse(File.open(File.join(TEST_DATA_DIR, 'minimal3.pep.xml')))
    pepxml.kind_of?(Bio::PepXML).should == true

    prot1 = pepxml.protein_name_to_object['>38SUR_2379_1524213_2']
    prot1.peptides.length.should == 2
    prot1.unique_spectra.should == 1
    prot1.non_unique_spectra.should == 1
    prot1.estimated_spectral_count.should == 2.0
  end

  it 'should parse when the protein and protein_desc attributes are both defined' do
    pepxml = Bio::PepXML.parse(File.open(File.join(TEST_DATA_DIR, 'contaminant.pep.xml')))
    pepxml.kind_of?(Bio::PepXML).should == true

    prot1 = pepxml.protein_name_to_object['CNTM:cont_sp']
    prot1.nil?.should == false
    prot1.identifier.should == 'CNTM:cont_sp'
    prot1.descriptive_name.should == 'CNTM:cont_sp P13647 K2C5_HUMAN Keratin, type II cytoskeletal 5 (Cytokeratin 5) (K5) (CK 5) (58 kDa cytokeratin) - Homo sapiens (Human). # pI:8.14 MW:62462'
  end
end
