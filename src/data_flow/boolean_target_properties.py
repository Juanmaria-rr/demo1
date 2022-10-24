from pyspark import F

### Still need the priotization of HPA main terms and predicted. 

def biotype_query(target, queryset):
    target_biotype = (target
    .select("id", "biotype")

    .dropDuplicates(['id','biotype'])
    .join(queryset, target.id == queryset.targetid, "right")
    .select("targetid", "biotype")
    .withColumn('Nr_biotype',
    F.when (F.col('biotype')=='protein_coding',1)
    .otherwise(0))
    )
    return target_biotype #### [id,biotype]

def target_membrane(target, queryset):
    
    membrane_terms=(parent_child_cousins_2
    .filter(F.col('Name') == 'Cell membrane')
    .select(parent_child_cousins_2.toSearch).rdd.flatMap(lambda x: x).collect()[0]
    )
    ###
    ###
    
    column= (target
    .select(
        F.col('id').alias('loc_id'),
        F.col('subcellularLocations'),
        F.explode_outer('subcellularLocations.termSL').alias('termSL')
    )
    .withColumn('isinMb', F.col('termSL').isin(membrane_terms))
    .withColumn('IsinMembrane', F.when(F.col('isinMb') == 'true', 'Yes')
    .when(F.col('isinMb') == 'false', 'No'))
    .filter(F.col('IsinMembrane')=='Yes')
    .dropDuplicates(['loc_id','IsinMembrane'])
    .select('loc_id','IsinMembrane')
    .withColumn('Nr_mb',
    F.when (F.col('IsinMembrane')=='Yes',F.lit(1))
    .when(F.col('IsinMembrane') == 'No',F.lit(0))
    .otherwise(None))
    .join(queryset, F.col('loc_id') == queryset.targetid, "right")
    #### conversion Membrane = 1
  
    )
    return column ### AND JOIN TO COLUMN.

#### 
def drug_query(molecule,molecule_mec, queryset):
    drug_approved=(molecule
    .select(
        F.col('id').alias('drug_id'),
        'isApproved',
        'linkedTargets')
    )
    drug_action=(molecule_mec
    .select('actionType', F.explode_outer(F.col('chemblIds')).alias('chembl')))

    appdrug_targets=(drug_approved
    .join(drug_action, drug_action.chembl == drug_approved.drug_id, 'left')
    .withColumn('targets',F.explode_outer(F.col('linkedTargets.rows')))
    .filter(F.col('isApproved')=='true')
    .select('targets','chembl','actionType')

    .dropDuplicates(['targets','chembl'])
    .groupBy('targets')
    .agg(F.collect_list('chembl').alias('App_drug_ChEMBL'),F.collect_list('actionType').alias('App_drug_actionType'))

    ### We add a new columns with Yes because all of them have a drug
    .withColumn('Drug',F.lit('Yes'))
    .withColumn('Nr_drug',
    F.when(F.col('Drug') !="",F.lit(1))
    .otherwise(None))
    .join(queryset, queryset.targetid == F.col('targets'),'right')
    )
    ### We put Yes/No in a new column 'boolean_drug' using all
    ## Basicaly we put No in the null values 
    return appdrug_targets 

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
    ## Make a column only for tag 'Yes'
    .withColumn('Chprob', F.lit('Yes'))
    .withColumn('Nr_chprob',
    F.when(F.col('Chprob') !="",F.lit(1))
    .otherwise(None))
    .join(queryset, F.col('chemid') == queryset.targetid, 'right')
    ## Tag those with something different to '' as No.
    )
    return chprob


#####

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
        'selection', F.when(F.col('cal_score') < (-0.1), str('-1'))
    .when(F.col('cal_score').between(-0.0999999,0.0999999),str('0'))
    .when(F.col('cal_score') > (0.1),str('1')))
    ## JOIN with Queryset
    .select('id','cal_score', 'selection', 'constraintType')
    .join(queryset, queryset.targetid == F.col('id'), 'right')
    )
    return loftolerance


##### 

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
    .withColumn('Nr_Mousemodels',
    F.when(F.col('Nr_mouse_models') !="0",F.lit(1))
    .otherwise(None))    
    .join(chemi_probes, chemi_probes.targetid == F.col('target_id_'),'right')

    )
    return moclass