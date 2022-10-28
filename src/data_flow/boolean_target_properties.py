from pyspark import F


def biotype_query(target, queryset):
    target_biotype = (
        target.select("id", "biotype")
        .dropDuplicates(["id", "biotype"])
        .join(queryset, target.id == queryset.targetid, "right")
        .select("targetid", "biotype")
        .withColumn(
            "Nr_biotype", F.when(F.col("biotype") == "protein_coding", 1).otherwise(0)
        )
    )
    return target_biotype


def target_membrane(target, queryset):  ### to solve 0, nulls and 1

    membrane_terms = (
        parent_child_cousins.filter(F.col("Name") == "Cell membrane")
        .select(parent_child_cousins.toSearch)
        .rdd.flatMap(lambda x: x)
        .collect()[0]
    )
    source_list = ["HPA_1", "HPA_secreted", "HPA_add_1", "uniprot_1"]
    membrane = (
        target.select(
            F.col("id").alias("loc_id"), F.explode_outer("subcellularLocations")
        )
        .select("loc_id", "col.*")
        .withColumn(
            "Nr_mb",
            F.when(
                (F.col("source") == "HPA_main")
                & (F.col("termSL").isin(membrane_terms)),
                F.lit("HPA_1"),
            )
            .when(
                (F.col("source") == "HPA_extracellular_location"), F.lit("HPA_secreted")
            )
            .when(
                (F.col("source") == "HPA_additional")
                & (F.col("termSL").isin(membrane_terms)),
                F.lit("HPA_add_1"),
            )
            .when(
                (F.col("source") == "uniprot") & (F.col("termSL").isin(membrane_terms)),
                F.lit("uniprot_1"),
            )
            .otherwise(F.lit(None)),
        )
        .filter(F.col("Nr_mb").isin(source_list))
        .select("loc_id", "Nr_mb", "source")
        .dropDuplicates(["loc_id", "Nr_mb"])
        .groupBy("loc_id")
        .agg(
            F.array_distinct(F.collect_list("Nr_mb")).alias("mb"),
            F.count("source").alias("counted"),
        )
        .withColumn(
            "loc",
            F.when(
                (F.array_contains(F.col("mb"), "HPA_secreted"))
                & (F.col("counted") == 1),
                "onlySecreted",
            )
            .when(
                (F.array_contains(F.col("mb"), "HPA_secreted"))
                & (F.col("counted") != 1),
                "secreted&inMembrane",
            )
            .otherwise(F.lit("inMembrane")),
        )
        .join(queryset, F.col("loc_id") == queryset.targetid, "right")
        .withColumn(
            "Nr_mb",
            F.when((F.col("loc") == "secreted&inMembrane"), F.lit(1))
            .when((F.col("loc") == "inMembrane"), F.lit(1))
            .otherwise(F.lit("0")),
        )
        .withColumn(
            "Nr_secreted",
            F.when((F.col("loc") == "secreted&inMembrane"), F.lit(1))
            .when((F.col("loc") == "onlySecreted"), F.lit(1))
            .otherwise(F.lit("0")),
        )
    )
    return membrane


def ligand_pocket_query(target, queryset):
    ligpock = (
        target.select(
            F.col("id").alias("target_id"),
            F.explode_outer("tractability").alias("new_struct"),
        )
        .filter(
            (
                (F.col("new_struct.id") == "High-Quality Ligand")
                & (F.col("new_struct.value") == True)
            )
            | (
                (F.col("new_struct.id") == "High-Quality Pocket")
                & (F.col("new_struct.value") == True)
            )
        )
        .withColumn("type", F.col("new_struct").getItem("id"))
        .join(queryset, F.col("target_id") == queryset.targetid, "right")
        #### .select("target_id", 'targetid', "type")
        .withColumn(
            "Nr_Ligand",
            F.when((F.col("type")) != "High-Quality Ligand", F.lit(1)).otherwise(
                F.lit(0)
            ),
        )
        .withColumn(
            "Nr_Pocket",
            F.when((F.col("type")) != "High-Quality Pocket", F.lit(1)).otherwise(
                F.lit(0)
            ),
        )
    )
    return ligpock


############


def safety_query(target, queryset):
    safety = (
        target.select(F.col("id").alias("saf_id"), F.explode_outer("safetyLiabilities"))
        .select("saf_id", F.col("col.*"))
        .groupBy("saf_id")
        .agg(
            F.count(F.col("event")).alias("nEvents"),
            F.array_distinct(F.collect_list("event")).alias("events"),
        )
        .withColumn(
            "nrEvent", F.when(F.col("nEvents") != 0, F.lit(1)).otherwise(F.lit(0))
        )
        ### Make the join
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


def paralogs(target, queryset):

    paralog = (
        target.select(F.col("id"), F.explode(F.col("homologues")))
        .select(F.col("id").alias("paralog_id"), F.col("col.*"))
        .withColumn("homoType", F.split(F.col("homologyType"), "_").getItem(0))
        .withColumn("howmany", F.split(F.col("homologyType"), "_").getItem(1))
        .withColumn("homoType", F.regexp_replace("homoType", "other", "paralog_other"))
        .withColumn(
            "homoType", F.regexp_replace("homoType", "within", "paralog_intrasp")
        )
        .select("paralog_id", "homologyType", "homoType", "howmany")
        .filter(F.col("homoType").contains("paralog"))
        .groupBy("paralog_id")
        .pivot("homoType")
        .agg(F.count("homoType"))
        .join(queryset, queryset.targetid == F.col("paralog_id"), "right")
        .withColumn(
            "nr_paralogs",
            F.when(
                (F.col("paralog_intrasp") > 0) | (F.col("paralog_other") > 0), F.lit(1)
            ).otherwise(F.lit(0)),
        )
    )
    return paralog


def orthologs_mouse(target, queryset):
    ### ONLY MOUSE orthologs
    ortholog = (
        target.select(F.col("id"), F.explode(F.col("homologues")))
        .select(F.col("id").alias("ortholog_id"), F.col("col.*"))
        .withColumn("homoType", F.split(F.col("homologyType"), "_").getItem(0))
        .withColumn("howmany", F.split(F.col("homologyType"), "_").getItem(1))
        .withColumn("homoType", F.regexp_replace("homoType", "other", "paralog_other"))
        .withColumn(
            "homoType", F.regexp_replace("homoType", "within", "paralog_intrasp")
        )
        .select(
            "ortholog_id",
            "homologyType",
            "homoType",
            "howmany",
            "targetGeneid",
            "targetPercentageIdentity",
            "queryPercentageIdentity",
            "speciesName",
        )
        .filter(
            (F.col("homoType").contains("ortholog")) & (F.col("speciesName") == "Mouse")
        )
        .groupBy("ortholog_id", "howmany")
        .pivot("homoType")
        .agg(
            F.expr("percentile_approx(targetPercentageIdentity, 0.5)").alias(
                "Target_median"
            ),
            F.expr("percentile_approx(queryPercentageIdentity, 0.5)").alias(
                "Query_median"
            ),
        )
        .withColumn("Nr_orth_Target_median", F.col("ortholog_Target_median") / 100)
        .withColumn("Nr_orth_Query_median", F.col("ortholog_Query_median") / 100)
        .join(queryset, queryset.targetid == F.col("ortholog_id"), "right")
    )
    return ortholog


def driver_genes(target, queryset):

    oncotsg_list = [
        "TSG",
        "oncogene",
        "Oncogene",
        "oncogene",
        "oncogene,TSG",
        "TSG,oncogene",
        "fusion,oncogene",
        "oncogene,fusion",
    ]

    driver = (
        target.select(
            "id", "approvedSymbol", F.explode_outer(F.col("hallmarks.attributes"))
        )
        .select(F.col("id").alias("driver_id"), F.col("col.description"))
        .withColumn(
            "Nr_CDG",
            F.when(F.col("description").isin(oncotsg_list), F.lit(1)).otherwise(
                F.lit(0)
            ),
        )
        .join(queryset, queryset.targetid == F.col("driver_id"), "right")
    )
    return driver


#############


def tep_query(target, queryset):
    ### Has TEP
    tep = (
        target.select(F.col("id").alias("tep_id"), F.col("tep.*"))
        .withColumn(
            "hasTEP", F.when(F.col("description") != "null", "Yes").otherwise("No")
        )
        .withColumn(
            "Nr_TEP",
            F.when(F.col("description") != "null", F.lit(1)).otherwise(F.lit(0)),
        )
        ### Make join
        .join(queryset, queryset.targetid == F.col("tep_id"), "right")
    )
    return tep


def mousemod_class(mouse, queryset):
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
        .join(queryset, queryset.targetid == F.col("target_id_"), "right")
        .withColumn(
            "Nr_Mousemodels",
            F.when(F.col("Nr_mouse_models") != "0", F.lit(1)).otherwise(0),
        )
    )
    return moclass


def chemical_probes(target, queryset):
    chprob = (
        target.select(F.col("id").alias("chemid"), F.explode(F.col("chemicalProbes")))
        .filter(F.col("col.isHighQuality") == "true")
        .select(F.col("chemid"), "col.*")
        .select(F.col("*"), F.explode(F.col("urls")))
        .select(F.col("*"), "col.*")
        .select("chemid", "mechanismOfAction")
        .groupBy("chemid", "mechanismOfAction")
        .agg(F.count("mechanismOfAction").alias("counts"))
        .select(
            "chemid",
            "counts",
            F.concat_ws(":", F.col("mechanismOfAction"), F.col("counts")).alias(
                "counted"
            ),
        )
        .where("counted!='0' ")
        .groupBy("chemid")
        .agg(
            F.collect_list("counted").alias("ChemicalProbes_HC"),
            F.count("counts").alias("count_chprob"),
        )
        .join(queryset, queryset.targetid == F.col("chemid"), "right")
        .withColumn(
            "Nr_chprob",
            F.when(F.col("count_chprob") != 0, F.lit(1)).otherwise(F.lit(0)),
        )
    )
    return chprob


def clin_trials(molecule, molecule_mec, queryset):
    clin_trials = [0, 1, 2, 3]
    drug_approved = molecule.select(
        F.col("id").alias("drug_id"), "maximumClinicalTrialPhase", "linkedTargets"
    ).withColumn(
        "Nr_ClinTrial",
        F.when(F.col("maximumClinicalTrialPhase") == 4, F.lit(1))
        .when(F.col("maximumClinicalTrialPhase").isin(clin_trials), F.lit(0))
        .otherwise(F.lit(None)),
    )
    drug_action = molecule_mec.select(
        "actionType", F.explode_outer(F.col("chemblIds")).alias("chembl")
    )

    appdrug_targets = (
        drug_approved.join(
            drug_action, drug_action.chembl == drug_approved.drug_id, "left"
        )
        .withColumn(
            "ClinTrials",
            F.concat_ws(
                "_",
                F.col("chembl"),
                F.lit("ClinTrialPhase"),
                F.col("maximumClinicalTrialPhase"),
            ),
        )
        .withColumn("targets", F.explode_outer(F.col("linkedTargets.rows")))
        .select(
            "targets",
            "chembl",
            "actionType",
            "maximumClinicalTrialPhase",
            "ClinTrials",
            "Nr_ClinTrial",
        )
        .dropDuplicates(["targets", "chembl"])
        .groupBy("targets")
        .agg(
            F.collect_list("ClinTrials").alias("ChEMBL&ClinTrialPhase"),
            F.collect_list("actionType").alias("App_drug_actionType"),
            F.max("Nr_ClinTrial").alias("maxClinTrialPhase"),
        )
        .join(queryset, queryset.targetid == F.col("targets"), "right")
    )
    ## Basicaly we put No in the null values
    return appdrug_targets
