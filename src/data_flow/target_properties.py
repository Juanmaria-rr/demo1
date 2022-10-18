from pyspark import F


def biotype_query(target, queryset):
    target_biotype = (target
    .select("id", "biotype")
    .join(queryset, target.id == queryset.targetid, "right"
    )
    .select("targetid", "biotype")
    )
    return target_biotype #### [id,biotype]

def target_location (target, queryset):
    column= (target
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

def drug_query(molecule,molecule_mec, queryset):
    drug_approved=(molecule
    .select('id','isApproved','linkedTargets')
    )
    drug_action=(molecule_mec
    .select('actionType', F.explode_outer(F.col('chemblIds')).alias('chembl')))

    appdrug_targets=(drug_approved
    .join(drug_action, drug_action.chembl == drug_approved.id, 'left')
    .withColumn('targets',F.explode_outer(F.col('linkedTargets.rows')))
    .filter(F.col('isApproved')=='true')
    .select('targets','chembl','actionType')
    .dropDuplicates(['targets','chembl'])
    .groupBy('targets')
    .agg(F.collect_list('chembl'),F.collect_list('actionType'))
    .join(queryset, queryset.targetid == F.col('targets'),'right')
    )
    return appdrug_targets 

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

def constraint (target,queryset):

    loftolerance=(target
    .select(
        F.col('id'),
        F.explode('constraint'))
    .select(
        F.col('id'),
        F.col('col.*'))
    .filter(
        F.col('constraintType')== 'lof')
    .withColumn(
        "cal_score", (F.col("upperRank") - 9456) / 19196)
    .withColumn(
        'selection', F.when(F.col('cal_score') < (-0.1), str('lofIntolerant')))
    .withColumn(
        'selection', F.when(F.col('cal_score').between(-1,-0.1),str('LOFintolerant'))
    .when(F.col('cal_score').between(-0.0999999,0.0999999),str('Neutral'))
    .when(F.col('cal_score').between(0.1,1),str('LOFtolerant')))
## falta hacer el join con el queryset de targets para darle el score de seleccion  
    .select('id','cal_score', 'selection', 'constraintType')
    .join(queryset, queryset.targetid == F.col('id'), 'right')
    )
    return loftolerance

    