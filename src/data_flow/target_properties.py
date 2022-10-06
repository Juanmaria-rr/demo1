def biotype_query(target, queryset):
    target_biotype = (target
    .select("id", "biotype")
    .join(queryset, target.id == queryset.targetid, "right"
    )
    .select("targetid", "biotype")
    )
    return target_biotype #### [id,biotype]

def target_location (target, queryset):
    column= (
    target
    .select('id',
        F.col('subcellularLocations'),
        F.explode_outer('subcellularLocations')
    )
    .select('id',
        F.col('col.location')
    )
    .groupBy('id')
    .agg(F.collect_list('location').alias('location'))
    .join(queryset,queryset.targetid == F.col('id'),'right')
    )
    return column ### AND JOIN TO COLUMN.

def drug_query(target, queryset):
    tractab = (
    target
    .select(
        F.col("id").alias("target_id"),
        F.explode_outer("tractability").alias("new_struct")
        )
    .select("target_id", F.col("new_struct.*"))
    .filter(
        (F.col("id") == "Approved Drug")
         & (F.col("value") == "True")
    )
    .groupBy('target_id').agg(F.collect_list('modality').alias('Approved_drugType'))
    .join(queryset, queryset.targetid == F.col('target_id'),'right')
    )
    return tractab 

def partner_drugs (molecule,interact_db,queryset): 
    tar_group=(molecule
    .select(
        F.col('id').alias('CHEMBL'), 
        F.explode(F.col('linkedTargets.rows')).alias('target'))
    .groupBy('target')
    .agg(
        F.collect_list('CHEMBL').alias('CHEMBL_list'),
        F.count('CHEMBL').alias('N_CHEMBL'))                                            
    )

    cnt_cond = lambda cond: F.sum(F.when(cond, 1).otherwise(0))
    ## Function to convert any number under condition to 1
    partner_drugs=(interact_db
    .filter(
        F.col('sourceDatabase') =='intact')
    .select(
        'sourceDatabase', 
        'targetA',
        'targetB',
        'scoring')
    .filter(
        F.col('scoring') > '0.42')
    .dropDuplicates(
        ["targetA","targetB"])
    .join(tar_group ,F.col('targetB') == tar_group.target,"left")
    .groupBy('targetA')
    .agg(
        cnt_cond(F.col('N_CHEMBL') > 0).alias('N_partner_drug'))
    .join(queryset ,queryset.targetid ==  F.col('targetA'),"right")    
    )
    return partner_drugs

def chemical_probes (target,queryset): 
    chprob=(target
    .select(
        F.col('id').alias('chemid'),
        F.explode(F.col('chemicalProbes')))
    .filter(F.col('col.isHighQuality') == 'true')
    .select(F.col('chemid'),'col.*')
    .select(
        F.col('*'),
        F.explode(F.col('urls')))
    .select(F.col('*'),'col.*')
    .select('chemid','mechanismOfAction')
    .groupBy('chemid','mechanismOfAction')
    .agg(
        F.count('mechanismOfAction').alias('counts'))
    .select(
        'chemid', 
        F.concat_ws(':',F.col('mechanismOfAction'),F.col('counts')).alias('counted'))
    .where("counted!='0' ")
    .groupBy('chemid')
    .agg(
        F.collect_list('counted').alias('ChemicalProbes_HC'))
    .join(queryset, F.col('chemid') == queryset.targetid, 'right')
    ) 
    return chprob 

def mousemod_class (mouse,chemi_probes): 

    moclass=(mouse
    .select(
        F.col('targetFromSourceId'),
        F.explode(F.col('modelPhenotypeClasses')).alias('classes'),
        F.col('classes.label'))
    .select(
        F.col('targetFromSourceId').alias('target_id_'), 
        F.col('classes.label'))
    .groupBy('target_id_')
    .agg(
        F.count('label').alias('Nr_mouse_models'),
        F.collect_set('label').alias('Different_PhenoClasses'))
    .join(chemi_probes, chemi_probes.targetid == F.col('target_id_'),'right')
    )
    return moclass 