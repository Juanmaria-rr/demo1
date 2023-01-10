from pyspark import F


def biotype_query(target, queryset):
    target_biotype = (
        target.select(
            F.col("id").alias('biotypeid'), "biotype")
        .dropDuplicates(["biotypeid", "biotype"])
        .select("biotypeid", "biotype")
        .withColumn('isProteinCoding', 
            F.when(F.col('biotype')=='protein_coding', F.lit('Yes'))
            .otherwise(F.lit('No')))
        .withColumn(
            "Nr_biotype", 
            F.when(F.col("biotype") == "protein_coding", F.lit(1))
            .otherwise(0))
        .join(queryset, F.col('biotypeid') == queryset.targetid, "right")
    )
    return target_biotype


def target_membrane(target, queryset):  ### to solve 0, nulls and 1

### Fine tunned function for including the secreted terms from uniprot.
### 03.11.2022 ###
    membrane_terms = (parent_child_cousins
        .filter(F.col("Name") == "Cell membrane")
        .select(parent_child_cousins.toSearch)
        .rdd.flatMap(lambda x: x)
        .collect()[0])

    secreted_terms =(parent_child_cousins
        .filter(F.col("Name") == "Secreted")
        .select(parent_child_cousins.toSearch)
        .rdd.flatMap(lambda x: x)
        .collect()[0])

    location_info= (target
    .select(
        F.col("id").alias("location_id"), F.explode_outer("subcellularLocations"))
    .withColumn('result', ## annotate which rows are 'null' 
        F.when (F.col('col.location').isNull(),'noInfo')
        .otherwise('hasInfo'))
    .select('location_id','result')
    .dropDuplicates(["location_id",'result'])
    )

    source_list = ["HPA_1", "HPA_secreted", "HPA_add_1", "uniprot_1",'uniprot_secreted']

    membrane = (
    target.select(
        F.col("id").alias("loc_id"), F.explode_outer("subcellularLocations")
    )
    .select("loc_id", "col.*")
    .withColumn(
        "Count_mb",
        F.when(
            (F.col("source") == "HPA_main")
            & (F.col("termSL").isin(membrane_terms)),
            F.lit("HPA_1")
        )
        .when(
            (F.col("source") == "HPA_extracellular_location"), 
            F.lit("HPA_secreted")
        )
        .when(
            (F.col("source") == "HPA_additional")
            & (F.col("termSL").isin(membrane_terms)),
            F.lit("HPA_add_1")
        )
        .when(
            (F.col("source") == "uniprot") & 
            (F.col("termSL").isin(membrane_terms)),
            F.lit("uniprot_1")
        )
        .when( ### Hacer el nuevo subset para 
            (F.col("source") == "uniprot")
            & (F.col("termSL").isin(secreted_terms)),
            F.lit("uniprot_secreted")
        )
        .otherwise(F.lit('Noinfo')),
    )
    .filter(F.col("Count_mb").isin(source_list))
    .select("loc_id", "Count_mb", "source")
    .dropDuplicates(["loc_id", "Count_mb"])
    .groupBy("loc_id")
    .agg(
        F.array_distinct(F.collect_list("Count_mb")).alias("mb"),
        F.count("source").alias("counted")
    )
    .withColumn( "loc",
        F.when(
            ((F.array_contains(F.col("mb"), "HPA_secreted") 
            & F.array_contains(F.col("mb"), "uniprot_secreted")
            ))
            & (F.col("counted") == 2), F.lit("onlySecreted")
        )
        .when(
            ((F.array_contains(F.col("mb"), "HPA_secreted") 
            | F.array_contains(F.col("mb"), "uniprot_secreted")
            ))
            & (F.col("counted") == 1),F.lit( "onlySecreted")
        )
        .when(
            ((F.array_contains(F.col("mb"), "HPA_secreted") &
            F.array_contains(F.col("mb"), "uniprot_secreted")))
            & (F.col("counted") > 2), F.lit("secreted&inMembrane"))
        .otherwise(F.lit("inMembrane"))
    )
    .join(location_info, F.col('loc_id')==location_info.location_id, 'right')
    .withColumn(
        "Nr_mb",
        F.when(
            (
                F.col("loc") == "secreted&inMembrane") | 
            (
                F.col("loc") == "inMembrane"),
            F.lit(1))
        .when(
            (
                (F.col("loc") != "inMembrane") | 
                (F.col("loc") != "secreted&inMembrane")
            ) &
            (F.col('result') == 'hasInfo'), F.lit(0))
        .when(F.col("result") == "noInfo", F.lit(None))
        .otherwise(F.lit(0))
        )
    .withColumn(
        "Nr_secreted",
        F.when(
            (
                F.col("loc") == "secreted&inMembrane") | 
            (
                F.col("loc") == "onlySecreted"),
            F.lit(1))
        .when(
            (
                (F.col("loc") != "onlySecreted") | 
                (F.col("loc") != "secreted&inMembrane")
            ) &
            (F.col('result') == 'hasInfo'), F.lit(0))
        .when(F.col("result") == "noInfo", F.lit(None))
        .otherwise(F.lit(0))
        )
    .withColumn(
        "name_mb",
        F.when(
            F.col("Nr_mb") == 1,
            F.lit('Yes'))
        .when(
            F.col("Nr_mb") == 0, F.lit('No'))
        .when(F.col("Nr_mb") == None, F.lit(None))
        )
    .withColumn(
        "name_secreted",
        F.when(
            F.col("Nr_secreted") == 1,
            F.lit('Yes'))
        .when(
            F.col("Nr_secreted") == 0, F.lit('No'))
        .when(F.col("Nr_secreted") == None, F.lit(None))
        )
    .join(queryset, F.col("loc_id") == queryset.targetid, "right")
    )
    return membrane

######### ----- ########
### fixed 04.11.2022
def ligand_pocket_query(target, queryset): 
    ligpock=(target
        .select(
            F.col("id").alias("target_id"),
            F.explode_outer("tractability").alias("new_struct"),
        )
        .filter(
            (
                (F.col("new_struct.id") == "High-Quality Ligand")
            )                        | 
                (
                (F.col("new_struct.id") == "High-Quality Pocket")
            )
        )
        .withColumn("type", F.col("new_struct").getItem("id"))
        .withColumn("presence", F.col("new_struct").getItem("value").cast('integer')) ## cast convert True = 1/ False = 0.
        .groupBy('target_id')
        .pivot('type')
        .agg(F.sum('presence'))
        .withColumn("Nr_Ligand", 
            F.when (
                F.col("High-Quality Ligand")== 1, F.lit(1))
            .otherwise(F.lit(0)))
        .withColumn("Nr_Pocket", 
            F.when (
                F.col("High-Quality Pocket")== 1, F.lit(1))
            .otherwise(F.lit(0)))
        .join(queryset, F.col("target_id") == queryset.targetid, "right")
        .withColumn('hasLigand', ### F.col('Nr_Ligand').cast('boolean'), 
            F.when(F.col('Nr_Ligand')==1, F.lit('Yes'))
            .when(F.col('Nr_Ligand')==0,F.lit('No'))
            .otherwise(F.lit(None))).na.fill("noInfo")
        .withColumn('hasPocket', ### F.col('Nr_Ligand').cast('boolean'), 
            F.when(F.col('Nr_Pocket')==1, F.lit('Yes'))
            .when(F.col('Nr_Pocket')==0,F.lit('No'))
            .otherwise(F.lit(None))).na.fill("noInfo")
        )
    return ligpock

############

def safety_query(target, queryset): ### 04.11.2022 donde anadimos columna 'info' con noReported para tagets con Null
    safety = (target
        .withColumn('info',
        F.when (F.col('safetyLiabilities')!=F.array(), F.lit('conInfo'))
        .otherwise(F.lit('noReported')))
        .select(F.col("id").alias("saf_id"), F.explode_outer("safetyLiabilities"),'info')
        .groupBy("saf_id",'info')
        .agg(
            F.count(F.col("col.event")).alias("nEvents"),
            F.array_distinct(F.collect_list("col.event")).alias("events"),
        )
        .withColumn(
            "hasSafetyEvent",
            F.when(
                (F.col("nEvents") > 0) 
                & (F.col('info')=='conInfo'),
                F.lit('Yes'))               
            .otherwise(F.lit(None))
        )
        .withColumn(
            "Nr_Event", 
            F.when(
                F.col("hasSafetyEvent")=='Yes', F.lit(-1))  
            .otherwise(F.lit(None))
        )
        .join(queryset, F.col("saf_id") == queryset.targetid, "right")
    )
    return safety


def constraint(target, queryset):  ### pending: solving hard coding of statistics

    loftolerance = (
        target.select(F.col("id").alias("constr_id"), F.explode("constraint"))
        .select(F.col("constr_id"), F.col("col.*"))
        .filter(F.col("constraintType") == "lof")
        .withColumn("cal_score", (F.col("upperRank") - 9456) / 19196)
        ## JOIN with Queryset
        .select("constr_id", "cal_score", "constraintType")
        .join(queryset, queryset.targetid == F.col("constr_id"), "right")
    )
    return loftolerance


def paralogs(target, queryset): ### HECHA columna Yes/No y cambiado orden de withColumn (la ponemos ahora antes del join)
                                ### 04.11.2022 Rehacer columna de numero porque tener paralogo es . 
                                ### Added continuous value from 60 to 100 and 0 for values from 0 to 60.
                                ### Negative value because is SAFETY
    paralog = (target
        .withColumn('hasInfo', 
            F.when(F.col('homologues')!=F.array(), F.lit('hasInfo'))
            .otherwise('noInfo/null'))
        .select(F.col("id").alias('paralog_id'),F.col('hasInfo'), F.explode_outer(F.col("homologues")),'hasInfo')
        .withColumn("homoType", F.split(F.col("col.homologyType"), "_").getItem(0))
        .withColumn("howmany", F.split(F.col("col.homologyType"), "_").getItem(1))
        .withColumn("homoType", F.regexp_replace("homoType", "other", "paralog_other"))
        .withColumn(
            "homoType", F.regexp_replace("homoType", "within", "paralog_intrasp")
        )
        .select("paralog_id","homoType", "howmany",'hasInfo','col.queryPercentageIdentity')
        .filter(F.col("homoType").contains("paralog"))
        .groupBy("paralog_id")
        .agg(
            F.max("queryPercentageIdentity").alias('max')
            )
        .withColumn('Nr_paralogs', 
        F.when(F.col('max')< 60, F.lit(0))
        .when(F.col('max') >=60, F.lit(-(((F.col('max')-60)*100)/40))))
        .join(queryset, queryset.targetid == F.col("paralog_id"), "right")
    )
    return paralog


def orthologs_mouse(target, queryset): ### Added continuous value from 80 to 100 and 0 for values from 0 to 80.
                                        #### Doability 
    ortholog = (target
        .select(F.col("id"), F.explode(F.col("homologues")))
        .select(F.col("id").alias("ortholog_id"), F.col("col.*"))
        .withColumn(
            "homoType", F.split(F.col("homologyType"), "_").getItem(0))
        .withColumn(
            "howmany", F.split(F.col("homologyType"), "_").getItem(1))
        .filter(
            (F.col("homoType").contains("ortholog")) 
            & (F.col("speciesName") == "Mouse") ### Filter by one_2_one? 
        )
        .select(
            "ortholog_id",
            "homoType", 
            "howmany",
            "targetGeneid",
            "targetPercentageIdentity",
            "queryPercentageIdentity"
        )
        .groupBy('ortholog_id')
        .agg(
            F.max('queryPercentageIdentity').alias('max'))
        .withColumn('Nr_ortholog', 
            F.when(F.col('max')< 80, F.lit(0))
            .when(F.col('max') >=80, F.lit((((F.col('max')-80)*100)/20))))  
        .join(queryset, queryset.targetid == F.col("ortholog_id"), "right")
        )
    return ortholog


def driver_genes(target, queryset): ## HECHA COLUMNA Yes/No
                                    ## Added -1/Null and Yes/Null 04.11.2022
    oncotsg_list = [                ## Fix duplicated targets 08.11.2022
        "TSG",
        "oncogene",
        "Oncogene",
        "oncogene",
        "oncogene,TSG",
        "TSG,oncogene",
        "fusion,oncogene",
        "oncogene,fusion",
    ]
    driver=(target
    .select(
        "id", "approvedSymbol", F.explode_outer(F.col("hallmarks.attributes")))
    .select(F.col("id").alias("driver_id"), F.col("col.description"))
    .withColumn('annotation', 
        F.when(F.col('description').isin(oncotsg_list), F.lit(1))
        .otherwise(F.lit(0)))
    .groupBy(F.col('driver_id'))
    .agg(
        F.max(F.col('annotation')).alias('counts')
    )
    .withColumn(
        "Nr_CDG",
        F.when(F.col("counts")!=0, F.lit(-1))
        .otherwise(
            F.lit(None)
        )
    )
    .withColumn(
        "isCancerDriverGene",
        F.when(F.col("counts")!=0, F.lit('Yes'))
        .otherwise(
            F.lit(None)
        )
    )
    .join(queryset, queryset.targetid == F.col("driver_id"), "right")
    )
    return driver


#############


def tep_query(target, queryset): ### Fixed with Null when not having TEP 04.11.2022
    ### Has TEP
    tep = (
        target.select(F.col("id").alias("tep_id"), F.col("tep.*"))
        .withColumn(
            "hasTEP", F.when(F.col("description") != "null", "Yes").otherwise("No")
        )
        .withColumn(
            "Nr_TEP",
            F.when(F.col("description") != "null", F.lit(1)).otherwise(F.lit(None)),
        )
        ### Make join
        .join(queryset, queryset.targetid == F.col("tep_id"), "right")
    )
    return tep


def mousemod_class(mouse, queryset): ### cambiado el orden de withColumn y agregada Yes/No
    moclass = (
        mouse.select(
            F.col("targetFromSourceId"),
            F.explode(F.col("modelPhenotypeClasses")).alias("classes"),
            F.col("classes.label"),
        )
        .select(F.col("targetFromSourceId").alias("target_id_"), F.col("classes.label"))
        .groupBy("target_id_")
        .agg(
            F.count("label").alias("Nr_mouse_models"),
            F.collect_set("label").alias("Different_PhenoClasses"),
        )
        .withColumn(
            "hasMouseKO",
            F.when(
                F.col("Nr_mouse_models") != "0", F.lit('Yes'))
                .otherwise(F.lit(0)),
        )        
        .withColumn(
            "Nr_Mousemodels",
            F.when(
                F.col("Nr_mouse_models") != "0", F.lit(1))
                .otherwise(F.lit(0)),
        )
        .join(queryset, queryset.targetid == F.col("target_id_"), "right")
    )
    return moclass


def chemical_probes(target, queryset): ### High-Quality Chemical Probes 
    chprob = (target                    ### Fixed 04.11.2022: No = does not have HQCP, Yes = it has, Null has no info.
        .select(
            F.col("id").alias("chemid"), 
            F.col("chemicalProbes"))
        .withColumn('info',
            F.when(F.col('chemicalProbes')!=F.array(),F.lit('hasInfo'))
            .otherwise(F.lit('noInfo')))
        .select(F.col("chemid"), F.explode_outer(F.col("chemicalProbes")),F.col('info'))
        .withColumn('Nr_chprob',
            F.when(
                (F.col('info')=='hasInfo')
                & (F.col('col.isHighQuality')=='true'),
                F.lit(1))
            .when(
                (F.col('info')=='hasInfo')
                & (F.col('col.isHighQuality')=='false'),
                F.lit(0))
            .otherwise(F.lit(None)))
        .groupBy('chemid')
        .agg(F.max(F.col('Nr_chprob')).alias('Nr_chprob'))
        .withColumn('hasHQCP',
            F.when(
                F.col('Nr_chprob')==1, F.lit('Yes'))
            .when(
                F.col('Nr_chprob')==0, F.lit('No'))
            .otherwise(F.lit(None)))
        #### .select('chemid','info','hasHQCP','Nr_HQCP')
        .join(queryset, queryset.targetid == F.col("chemid"), "right")
        )
    return chprob


def clin_trials(molecule, molecule_mec, queryset):
    clin_trials = [0, 1, 2, 3]
    drug_approved = (molecule ### Fixed. Number of Max Clin Trial Phase (from 0 to 4) and Null when not having. 04.11.2022
    .select(
        F.col("id").alias("drug_id"), "maximumClinicalTrialPhase", "linkedTargets")
    .withColumn("targets", F.explode_outer(F.col("linkedTargets.rows")))
    .select(
        "drug_id",
        "targets",
        "maximumClinicalTrialPhase"
    )
    .dropDuplicates(["targets", "maximumClinicalTrialPhase"])
    .groupBy("targets")
    .agg(
        F.max("maximumClinicalTrialPhase").alias("maxClinTrialPhase"),
    )
    .withColumn(
        "inClinicalTrials",
        F.when(F.col("maxClinTrialPhase").isin([0,1,2,3,4]), F.lit('inClinicalTrials'))
        .otherwise(F.lit(None)))
    ### .select('targets','maxClinTrialPhase')
    .join(queryset, queryset.targetid == F.col("targets"), "right")
    )
    return drug_approved


def tissue_specific (hpa_data, queryset): 

### On 01.11.2022 we add the tissue distribution column: 

    import json
    import pandas as pd
    from pyspark.sql import SQLContext     

    cols_of_interest = [
        "Ensembl",
        "RNA tissue distribution",
        "RNA tissue specificity",
        "Antibody"    
        ]

    ## hpa_data = 'proteinatlas.json'

    with open(hpa_data, 'r') as f:
        for line in f:
            entry = json.loads(line)
            df = (
                    pd.DataFrame(entry)
                    .filter(items=cols_of_interest)
            )
    hpa_df=spark.createDataFrame(df)

    hpa=(hpa_df
    .select('Ensembl','RNA tissue distribution','RNA tissue specificity', 'Antibody')
    .withColumnRenamed('RNA tissue distribution','Tissue_distribution_RNA')
    .withColumnRenamed('RNA tissue specificity','Tissue_specificity_RNA')
    .select('Ensembl','Tissue_specificity_RNA','Tissue_distribution_RNA')
    .withColumn('Nr_specificity', ### Fixed numbers asigned to the distribution thing (is SAFETY) 04.11.2022
        F.when(F.col('Tissue_specificity_RNA')=='Tissue enriched',F.lit(1))### tissue enriched 1
        .when(F.col('Tissue_specificity_RNA')=='Group enriched', F.lit(0.5)) ### group enriched 0.5 
        .when(F.col('Tissue_specificity_RNA')=='Tissue enhanced', F.lit(0.5))
        .when(F.col('Tissue_specificity_RNA')=='Low tissue specificity', F.lit(-1))
        .when(F.col('Tissue_specificity_RNA')=='Not detected', F.lit(None)))
    .withColumn('Nr_distribution', ### Fixed numbers asigned to the distribution thing (is SAFETY) 04.11.2022
    F.when(F.col('Tissue_distribution_RNA')=='Detected in single', F.lit(1)) ### single 1 // only in one tissue
    .when(F.col('Tissue_distribution_RNA')=='Detected in some', F.lit(0.5)) ### some 0.5 // more thann one but less than 1/3
    .when(F.col('Tissue_distribution_RNA')=='Detected in many', F.lit(0)) ## at least 1/3 but not all tissues
    .when(F.col('Tissue_distribution_RNA')=='Detected in all', F.lit(-1))
    .when(F.col('Tissue_distribution_RNA')=='Not detected', F.lit(None)))
    .join(queryset, queryset.targetid == F.col("Ensembl"), "right")
    )
    return hpa