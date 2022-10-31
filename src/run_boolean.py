from target_engine_repo.src.data_flow.boolean_target_properties import biotype_query, chemical_probes, clin_trials, constraint, driver_genes, ligand_pocket_query, mousemod_class, orthologs_mouse, paralogs, safety_query, target_membrane, tep_query, tissue_specific


queryset=target.select('id').withColumnRenamed('id','targetid').limit(61524)
hpa_data = 'proteinatlas.json'

biotype = biotype_query(target, queryset)
location = target_membrane(target, biotype)
ligand_pocket = ligand_pocket_query (target, location)
safety = safety_query (target, ligand_pocket)
selection = constraint(target, safety)
paralog = paralogs(target, selection)
ortholog = orthologs_mouse(target, paralog)
drivers = driver_genes(target, ortholog)
tep = tep_query(target, drivers)
mice = mousemod_class(mouse, tep)
chemprob = chemical_probes(target, mice)
drug_clin = clin_trials(molecule, molecule_mec, chemprob)
tissue_specificity = tissue_specific (hpa_data, chemprob)

#Selection of relevant columns
info=(drug_clin
    .select(
        'targetid', 
        'Nr_biotype',
        'Nr_mb',
        'Nr_secreted',
        'Nr_Pocket',
        'Nr_Ligand',
        'nrEvent',
        'cal_score',
        'nr_paralogs',
        'Nr_orth_Query_median',
        'Nr_CDG',
        'Nr_TEP',
        'Nr_Mousemodels',
        'Nr_chprob',
        'maxClinTrialPhase'
    )
    .withColumnRenamed('Nr_biotype', 'isProteinCoding')
    .withColumnRenamed('Nr_mb', 'isInMembrane')
    .withColumnRenamed('Nr_secreted', 'isSecreted')
    .withColumnRenamed('Nr_Pocket', 'hasPocket')
    .withColumnRenamed('Nr_Ligand', 'hasLigand')
    .withColumnRenamed('nrEvent', 'hasSafetyEvent')
    .withColumnRenamed('cal_score','geneticConstraint')
    .withColumnRenamed('nr_paralogs', 'hasParalogs')
    .withColumnRenamed('Nr_orth_Query_median','mouseOrthologIdentityPercentage')
    .withColumnRenamed('Nr_CDG','isCancerDriverGene')
    .withColumnRenamed('Nr_TEP','hasTEP')
    .withColumnRenamed('Nr_Mousemodels','hasMouseKO')
    .withColumnRenamed('Nr_chprob','hasHighQualityChemicalProbes')
    .withColumnRenamed('maxClinTrialPhase','inClinicalTrials')
)



