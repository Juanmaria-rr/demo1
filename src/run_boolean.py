from target_engine_repo.src.data_flow.boolean_target_properties import biotype_query
from target_engine_repo.src.data_flow.target_properties import drug_query, target_membrane


biotype = biotype_query(target, queryset)
location = target_membrane(target,biotype) 
drug = drug_query(molecule,molecule_mec, location)
chemi_probes= chemical_probes (target,drug)
mouse_models= mousemod_class (mouse,chemi_probes)
target_constraint = constraint (target,mouse_models)

#Selection of relevant columns
info=(target_constraint
    .select(
        'targetid', 
        'Nr_biotype',
        'Nr_mb',
        'Nr_drug', ## Any type of drug 'Approved' for the target
        'Nr_chprob', ## Type and number of High Confidence chemical probes 
        'Nr_Mousemodels', ## Number of mice models using containing the target
        ## Distinct classes of mice models phenotypes
        'cal_score' ##Classification of constraint towards Loss of function 
    )
    .withColumnRenamed('Nr_biotype', 'isProteinCoding')
    .withColumnRenamed('Nr_mb', 'isInMembrane')
    .withColumnRenamed('Nr_drug','hasDrugApproved')
    .withColumnRenamed('Nr_chprob','hasHighQualityChemicalProbes')
    .withColumnRenamed('Nr_Mousemodels','hasMouseKO')
    .withColumnRenamed('cal_score','geneticConstraint')
)