module Bio::DivvyProteomics::DivvyableProtein
  def unique_spectra
    return 0 if @peptides.nil? or @peptides.empty?
    num = @peptides.select{|pep| pep.parent_proteins.length == 1}.collect{|pep| pep.redundancy}.reduce(:+)
    num ||= 0
    return num
  end

  def non_unique_spectra
    return 0 if @peptides.nil? or @peptides.empty?
    num = @peptides.reject{|pep| pep.parent_proteins.length == 1}.collect{|pep| pep.redundancy}.reduce(:+)
    num ||= 0

    return num
  end

  # Are there any peptides that are assigned exclusively to this protein?
  def uniquely_identified_by_any_peptides?
    unique_spectra > 0
  end

  def estimated_spectral_count
    # How many unique spectra are there for each protein that shares a peptide with the current peptide
    return 0 if @peptides.nil? or @peptides.empty?
    peptide_shares = []
    # If all peptides are non-unique and shared with some number of other proteins, then output a negative number num shared spectra divided by the number of proteins
    if !uniquely_identified_by_any_peptides?
      # Don't attempt to divvy these up, because there are too many assumptions involved
      return 0
    else
      peptides.each do |peptide|
        log.debug "Tallying peptide #{peptide.identifier}, which is has #{peptide.redundancy} spectra shared among #{peptide.parent_proteins.length} proteins"
        log.debug "These proteins have #{peptide.parent_proteins.collect{|pro| pro.unique_spectra}.inspect} unique spectra each"
        total_linked_unique_spectra = peptide.parent_proteins.collect{|pro| pro.unique_spectra}.reduce(:+)
        peptide_shares.push unique_spectra.to_f/total_linked_unique_spectra*peptide.redundancy
      end
      return peptide_shares.reduce(:+)
    end
  end
end
