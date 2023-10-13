from pyspark.sql import Window
import pyspark.sql.functions as F
from pyspark.ml.functions import array_to_vector
from pyspark.sql.functions import udf


def biotype_query(target, queryset):
    target_biotype = (
        target.select(F.col("id").alias("biotypeid"), "biotype")
        .select("biotypeid", "biotype")
        .withColumn(
            "isProteinCoding",
            F.when(F.col("biotype") == "protein_coding", F.lit("Yes"))
            .when((F.col("biotype") == ""), F.lit("No_info"))
            .otherwise(F.lit("No")),
        )
        .withColumn(
            "Nr_biotype",
            F.when(F.col("biotype") == "protein_coding", F.lit(1))
            .when((F.col("biotype") == ""), F.lit(None))
            .otherwise(0),
        )
        .join(queryset, F.col("biotypeid") == queryset.targetid, "right")
    )
    return target_biotype


def target_membrane(target, queryset):

    membrane_terms = (
        parent_child_cousins.filter(F.col("Name") == "Cell membrane")
        .select(parent_child_cousins.toSearch)
        .rdd.flatMap(lambda x: x)
        .collect()[0]
    )

    secreted_terms = (
        parent_child_cousins.filter(F.col("Name") == "Secreted")
        .select(parent_child_cousins.toSearch)
        .rdd.flatMap(lambda x: x)
        .collect()[0]
    )

    location_info = (
        target.select(
            F.col("id").alias("location_id"), F.explode_outer("subcellularLocations")
        )
        .withColumn(
            "result",
            F.when(F.col("col.location").isNull(), "noInfo").otherwise("hasInfo"),
        )
        .select("location_id", "result")
        .dropDuplicates(["location_id", "result"])
    )

    source_list = [
        "HPA_1",
        "HPA_secreted",
        "HPA_add_1",
        "uniprot_1",
        "uniprot_secreted",
        "HPA_dif",
    ]
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
                F.lit("HPA_1"),
            )
            .when(
                (F.col("source") == "HPA_main")
                & (F.col("termSL").isin(membrane_terms) == False),
                F.lit("HPA_dif"),
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
                (F.col("source") == "HPA_additional")
                & (F.col("termSL").isin(membrane_terms) == False),
                F.lit("HPA_dif"),
            )
            .when(
                (F.col("source") == "uniprot") & (F.col("termSL").isin(membrane_terms)),
                F.lit("uniprot_1"),
            )
            .when(
                (F.col("source") == "uniprot") & (F.col("termSL").isin(secreted_terms)),
                F.lit("uniprot_secreted"),
            )
            .otherwise(F.lit("Noinfo")),
        )
        .filter(F.col("Count_mb").isin(source_list))
        .select("loc_id", "Count_mb", "source")
        .dropDuplicates(["loc_id", "Count_mb"])
        .groupBy("loc_id")
        .agg(
            F.collect_set("Count_mb").alias("mb"),
            F.count("source").alias("counted"),
        )
        .withColumn(
            "HPA_membrane",
            F.when(
                (
                    (F.array_contains(F.col("mb"), "HPA_1"))
                    | (F.array_contains(F.col("mb"), "HPA_add_1"))
                ),
                F.lit("yes"),
            )
            .when((F.array_contains(F.col("mb"), "HPA_dif")), F.lit("dif"))
            .otherwise(F.lit("no")),
        )
        .withColumn(
            "HPA_secreted",
            F.when(
                F.array_contains(F.col("mb"), "HPA_secreted"), F.lit("yes")
            ).otherwise(F.lit("no")),
        )
        .withColumn(
            "uniprot_membrane",
            F.when(F.array_contains(F.col("mb"), "uniprot_1"), F.lit("yes")).otherwise(
                F.lit("no")
            ),
        )
        .withColumn(
            "uniprot_secreted",
            F.when(
                F.array_contains(F.col("mb"), "uniprot_secreted"), F.lit("yes")
            ).otherwise(F.lit("no")),
        )
        .withColumn(
            "loc",
            F.when(
                ((F.col("HPA_membrane") == "yes")) & (F.col("HPA_secreted") == "no"),
                F.lit("inMembrane"),
            )
            .when(
                ((F.col("HPA_membrane") == "no") | (F.col("HPA_membrane") == "dif"))
                & (F.col("HPA_secreted") == "yes"),
                F.lit("onlySecreted"),
            )
            .when(
                (F.col("HPA_membrane") == "yes") & (F.col("HPA_secreted") == "yes"),
                F.lit("secreted&inMembrane"),
            )
            .when(
                (F.col("HPA_membrane") == "no") & (F.col("HPA_secreted") == "no"),
                F.when(
                    (F.col("uniprot_membrane") == "yes")
                    & (F.col("uniprot_secreted") == "no"),
                    F.lit("inMembrane"),
                )
                .when(
                    (F.col("uniprot_membrane") == "no")
                    & (F.col("uniprot_secreted") == "yes"),
                    F.lit("onlySecreted"),
                )
                .when(
                    (F.col("uniprot_membrane") == "yes")
                    & (F.col("uniprot_secreted") == "yes"),
                    F.lit("secreted&inMembrane"),
                ),
            )
            .when(F.col("HPA_membrane") == "dif", F.lit("noMembraneHPA")),
        )
        .join(queryset, F.col("loc_id") == queryset.targetid, "right")
        .join(location_info, F.col("targetid") == location_info.location_id, "left")
        .withColumn(
            "Nr_mb",
            F.when(
                (
                    (F.col("loc") == "secreted&inMembrane")
                    | (F.col("loc") == "inMembrane")
                ),
                F.lit(1),
            )
            .when(
                (
                    (F.col("loc") != "secreted&inMembrane")
                    | (F.col("loc") != "inMembrane")
                ),
                F.lit(0),
            )
            .when((F.col("loc").isNull()) & (F.col("result") == "hasInfo"), F.lit(0)),
        )
        .withColumn(
            "Nr_secreted",
            F.when(
                (F.col("loc") == "secreted&inMembrane")
                | (F.col("loc") == "onlySecreted"),
                F.lit(1),
            )
            .when(
                ((F.col("loc") == "inMembrane")) & (F.col("result") == "hasInfo"),
                F.lit(0),
            )
            .when(
                (
                    (F.col("loc") != "onlySecreted")
                    | (F.col("loc") != "secreted&inMembrane")
                )
                & (F.col("result") == "hasInfo"),
                F.lit(0),
            )
            .when(F.col("result") == "noInfo", F.lit(None))
            .otherwise(F.lit(0)),
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
            ((F.col("new_struct.id") == "High-Quality Ligand"))
            | ((F.col("new_struct.id") == "High-Quality Pocket"))
            | ((F.col("new_struct.id") == "Small Molecule Binder"))
        )
        .withColumn("type", F.col("new_struct").getItem("id"))
        .withColumn("presence", F.col("new_struct").getItem("value").cast("integer"))
        .groupBy("target_id")
        .pivot("type")
        .agg(F.sum("presence"))
        .withColumn(
            "Nr_Ligand",
            F.when(F.col("High-Quality Ligand") == 1, F.lit(1)).otherwise(F.lit(0)),
        )
        .withColumn(
            "Nr_Pocket",
            F.when(F.col("High-Quality Pocket") == 1, F.lit(1)).otherwise(F.lit(0)),
        )
        .withColumn(
            "Nr_sMBinder",
            F.when(F.col("Small Molecule Binder") == 1, F.lit(1)).otherwise(F.lit(0)),
        )
        .join(queryset, F.col("target_id") == queryset.targetid, "right")
        .withColumn(
            "hasLigand",
            F.when(F.col("Nr_Ligand") == 1, F.lit("Yes"))
            .when(F.col("Nr_Ligand") == 0, F.lit("No"))
            .otherwise(F.lit(None)),
        )
        .withColumn(
            "hasPocket",
            F.when(F.col("Nr_Pocket") == 1, F.lit("Yes"))
            .when(F.col("Nr_Pocket") == 0, F.lit("No"))
            .otherwise(F.lit(None)),
        )
    )
    return ligpock


def safety_query(target, queryset):
    safety = (
        target.withColumn(
            "info",
            F.when(F.col("safetyLiabilities") != F.array(), F.lit("conInfo")).otherwise(
                F.lit("noReported")
            ),
        )
        .select(
            F.col("id").alias("saf_id"), F.explode_outer("safetyLiabilities"), "info"
        )
        .groupBy("saf_id", "info")
        .agg(
            F.count(F.col("col.event")).alias("nEvents"),
            F.array_distinct(F.collect_list("col.event")).alias("events"),
        )
        .withColumn(
            "hasSafetyEvent",
            F.when(
                (F.col("nEvents") > 0) & (F.col("info") == "conInfo"), F.lit("Yes")
            ).otherwise(F.lit(None)),
        )
        .withColumn(
            "Nr_Event",
            F.when(F.col("hasSafetyEvent") == "Yes", F.lit(-1)).otherwise(F.lit(None)),
        )
        .join(queryset, F.col("saf_id") == queryset.targetid, "right")
    )
    return safety


def constraint(target, queryset):

    minUpperRank = (
        target.select(F.col("id").alias("constr_id"), F.explode("constraint"))
        .select(F.col("col.*"))
        .filter(F.col("constraintType") == "lof")
        .groupBy("constraintType")
        .agg(F.min("upperRank").alias("upperRank"))
        .select("upperRank")
        .rdd.flatMap(lambda x: x)
        .collect()[0]
    )

    maxUpperRank = (
        target.select(F.col("id").alias("constr_id"), F.explode("constraint"))
        .select(F.col("col.*"))
        .filter(F.col("constraintType") == "lof")
        .groupBy("constraintType")
        .agg(F.max("upperRank").alias("upperRank"))
        .select("upperRank")
        .rdd.flatMap(lambda x: x)
        .collect()[0]
    )

    loftolerance = (
        target.select(F.col("id").alias("constr_id"), F.explode("constraint"))
        .select(F.col("constr_id"), F.col("col.*"))
        .filter(F.col("constraintType") == "lof")
        .withColumn(
            "cal_score",
            F.lit(
                (
                    2
                    * (
                        (F.col("upperRank") - minUpperRank)
                        / (maxUpperRank - minUpperRank)
                    )
                )
                - 1
            ),
        )
        .select("constr_id", "cal_score", "constraintType")
        .join(queryset, queryset.targetid == F.col("constr_id"), "right")
    )
    return loftolerance


def paralogs(target, queryset):
    paralog = (
        target.withColumn(
            "hasInfo",
            F.when(F.col("homologues") != F.array(), F.lit("hasInfo")).otherwise(
                "noInfo/null"
            ),
        )
        .select(
            F.col("id").alias("paralog_id"),
            F.col("hasInfo"),
            F.explode_outer(F.col("homologues")),
            "hasInfo",
        )
        .withColumn("homoType", F.split(F.col("col.homologyType"), "_").getItem(0))
        .withColumn("howmany", F.split(F.col("col.homologyType"), "_").getItem(1))
        .withColumn("homoType", F.regexp_replace("homoType", "other", "paralog_other"))
        .withColumn(
            "homoType", F.regexp_replace("homoType", "within", "paralog_intrasp")
        )
        .select(
            "paralog_id",
            "homoType",
            "howmany",
            "hasInfo",
            "col.queryPercentageIdentity",
        )
        .filter(F.col("homoType").contains("paralog"))
        .groupBy("paralog_id")
        .agg(F.max("queryPercentageIdentity").alias("max"))
        .withColumn(
            "Nr_paralogs",
            F.when(F.col("max") < 60, F.lit(0)).when(
                F.col("max") >= 60, F.lit(-((F.col("max") - 60) / 40))
            ),
        )
        .join(queryset, queryset.targetid == F.col("paralog_id"), "right")
    )
    return paralog


def orthologs_mouse(target, queryset):
    ortholog = (
        target.select(F.col("id"), F.explode(F.col("homologues")))
        .select(F.col("id").alias("ortholog_id"), F.col("col.*"))
        .withColumn("homoType", F.split(F.col("homologyType"), "_").getItem(0))
        .withColumn("howmany", F.split(F.col("homologyType"), "_").getItem(1))
        .filter(
            (F.col("homoType").contains("ortholog")) & (F.col("speciesName") == "Mouse")
        )
        .select(
            "ortholog_id",
            "homoType",
            "howmany",
            "targetGeneid",
            "targetPercentageIdentity",
            "queryPercentageIdentity",
        )
        .groupBy("ortholog_id")
        .agg(F.max("queryPercentageIdentity").alias("max"))
        .withColumn(
            "Nr_ortholog",
            F.when(F.col("max") < 80, F.lit(0)).when(
                F.col("max") >= 80, F.lit((F.col("max") - 80) / 20)
            ),
        )
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
            "annotation",
            F.when(F.col("description").isin(oncotsg_list), F.lit(1)).otherwise(
                F.lit(0)
            ),
        )
        .groupBy(F.col("driver_id"))
        .agg(F.max(F.col("annotation")).alias("counts"))
        .withColumn(
            "Nr_CDG", F.when(F.col("counts") != 0, F.lit(-1)).otherwise(F.lit(None))
        )
        .withColumn(
            "isCancerDriverGene",
            F.when(F.col("counts") != 0, F.lit("Yes")).otherwise(F.lit(None)),
        )
        .join(queryset, queryset.targetid == F.col("driver_id"), "right")
    )
    return driver


def essentiality(geneEssentiality, queryset):
    essential = (
        geneEssentiality.select(
            "id", F.explode_outer("geneEssentiality").alias("geneEssentiality")
        )
        .withColumn(
            "Nr_essential", -F.col("geneEssentiality.isEssential").cast("integer")
        )
        .select(F.col("id").alias("idEssential"), "Nr_essential")
        .join(queryset, queryset.targetid == F.col("idEssential"), "right")
    )

    return essential


def tep_query(target, queryset):
    tep = (
        target.select(F.col("id").alias("tep_id"), F.col("tep.*"))
        .withColumn(
            "hasTEP", F.when(F.col("description") != "null", "Yes").otherwise("No")
        )
        .withColumn(
            "Nr_TEP",
            F.when(F.col("description") != "null", F.lit(1)).otherwise(F.lit(None)),
        )
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
        .withColumn(
            "hasMouseKO",
            F.when(F.col("Nr_mouse_models") != "0", F.lit("Yes")).otherwise(F.lit(0)),
        )
        .withColumn(
            "Nr_Mousemodels",
            F.when(F.col("Nr_mouse_models") != "0", F.lit(1)).otherwise(F.lit(0)),
        )
        .join(queryset, queryset.targetid == F.col("target_id_"), "right")
    )
    return moclass


def chemical_probes(target, queryset):
    chprob = (
        target.select(F.col("id").alias("chemid"), F.col("chemicalProbes"))
        .withColumn(
            "info",
            F.when(F.col("chemicalProbes") != F.array(), F.lit("hasInfo")).otherwise(
                F.lit("noInfo")
            ),
        )
        .select(
            F.col("chemid"), F.explode_outer(F.col("chemicalProbes")), F.col("info")
        )
        .withColumn(
            "Nr_chprob",
            F.when(
                (F.col("info") == "hasInfo") & (F.col("col.isHighQuality") == "true"),
                F.lit(1),
            )
            .when(
                (F.col("info") == "hasInfo") & (F.col("col.isHighQuality") == "false"),
                F.lit(0),
            )
            .otherwise(F.lit(None)),
        )
        .groupBy("chemid")
        .agg(F.max(F.col("Nr_chprob")).alias("Nr_chprob"))
        .withColumn(
            "hasHQCP",
            F.when(F.col("Nr_chprob") == 1, F.lit("Yes"))
            .when(F.col("Nr_chprob") == 0, F.lit("No"))
            .otherwise(F.lit(None)),
        )
        .join(queryset, queryset.targetid == F.col("chemid"), "right")
    )
    return chprob


def clin_trials(molecule, queryset):
    drug_approved = (
        molecule.select(
            F.col("id").alias("drug_id"), "maximumClinicalTrialPhase", "linkedTargets"
        )
        .withColumn("targets", F.explode_outer(F.col("linkedTargets.rows")))
        .select("drug_id", "targets", "maximumClinicalTrialPhase")
        .dropDuplicates(["targets", "maximumClinicalTrialPhase"])
        .groupBy("targets")
        .agg(
            F.max("maximumClinicalTrialPhase").alias("maxClinTrialPhase"),
        )
        .withColumn(
            "inClinicalTrials",
            F.when(
                F.col("maxClinTrialPhase") >= 0, F.lit(F.col("maxClinTrialPhase") / 4)
            ).otherwise(F.lit(None)),
        )
        .join(queryset, queryset.targetid == F.col("targets"), "right")
    )
    return drug_approved


def tissue_specific(hpa_data, queryset):

    cols_of_interest = [
        "Ensembl",
        "RNA tissue distribution",
        "RNA tissue specificity",
        "Antibody",
    ]

    ## hpa_data = 'proteinatlas.json'

    with open(hpa_data, "r") as f:
        for line in f:
            entry = json.loads(line)
            df = pd.DataFrame(entry).filter(items=cols_of_interest)
    hpa_df = spark.createDataFrame(df)

    hpa = (
        hpa_df.select(
            "Ensembl", "RNA tissue distribution", "RNA tissue specificity", "Antibody"
        )
        .withColumnRenamed("RNA tissue distribution", "Tissue_distribution_RNA")
        .withColumnRenamed("RNA tissue specificity", "Tissue_specificity_RNA")
        .select("Ensembl", "Tissue_specificity_RNA", "Tissue_distribution_RNA")
        .withColumn(
            "Nr_specificity",
            F.when(F.col("Tissue_specificity_RNA") == "Tissue enriched", F.lit(1))
            .when(F.col("Tissue_specificity_RNA") == "Group enriched", F.lit(0.75))
            .when(F.col("Tissue_specificity_RNA") == "Tissue enhanced", F.lit(0.5))
            .when(
                F.col("Tissue_specificity_RNA") == "Low tissue specificity", F.lit(-1)
            )
            .when(F.col("Tissue_specificity_RNA") == "Not detected", F.lit(None)),
        )
        .withColumn(
            "Nr_distribution",
            F.when(F.col("Tissue_distribution_RNA") == "Detected in single", F.lit(1))
            .when(F.col("Tissue_distribution_RNA") == "Detected in some", F.lit(0.5))
            .when(F.col("Tissue_distribution_RNA") == "Detected in many", F.lit(0))
            .when(F.col("Tissue_distribution_RNA") == "Detected in all", F.lit(-1))
            .when(F.col("Tissue_distribution_RNA") == "Not detected", F.lit(None)),
        )
        .join(queryset, queryset.targetid == F.col("Ensembl"), "right")
    )

    return hpa


def mouse_phenotypes(mousePhenotypes, queryset):

    mopheScore_path = "/Users/juanr/Desktop/target_engine_repo/src/data_flow/phenotypeScores/20230825_mousePheScores.csv"
    mopheScore = spark.read.csv(mopheScore_path, header=True)
    mousePhenoScoreFilter = mopheScore.select(
        F.col("id").alias("idLabel"),
        F.col("label").alias("phenoLabel"),
        F.col("score"),
    ).withColumn(
        "score",
        F.when(F.col("score") == 0.0, F.lit(0)).otherwise(F.lit(F.col("score"))),
    )

    ### Define Harmonic Sum function:
    def harmonic_sum(evidence_scores):
        harmonic_sum = sum(
            score / ((i + 1) ** (2)) for i, score in enumerate(evidence_scores)
        )
        return float(harmonic_sum)

    ### Define max Harmonic Sum function:
    def max_harmonic_sum(evidence_scores):
        max_theoretical_harmonic_sum = sum(
            1 / ((i + 1) ** (2)) for i in range(len(evidence_scores))
        )
        return float(max_theoretical_harmonic_sum)

    ### define function to scale the harmonic sum
    def scaledHarmonic(score, maximum):
        scaled_harmonic = score / maximum
        return float(scaled_harmonic)

    harmonic_sum_udf = udf(harmonic_sum)
    max_harmonic_sum_udf = udf(max_harmonic_sum)
    scaledHarmonic_udf = udf(scaledHarmonic)
    ### window function to take maximum of all harmonic sum
    window = Window.orderBy()

    scoreAggregation = (
        mousePhenotypes.select(
            "targetFromSourceId",
            F.explode_outer(F.col("modelPhenotypeClasses.id")).alias("id"),
        )
        .join(
            mousePhenoScoreFilter, F.col("id") == mousePhenoScoreFilter.idLabel, "left"
        )
        .withColumn("score", F.col("score").cast("float"))
        .groupBy("targetFromSourceId")
        .agg(array_to_vector(F.collect_list("score")).alias("score"))
        .withColumn("harmonic_sum", harmonic_sum_udf("score").cast("float"))
        .withColumn("maxHarmonicSum", max_harmonic_sum_udf("score").cast("float"))
        .withColumn("maximum", F.max("maxHarmonicSum").over(window).cast("float"))
        .withColumn("scaledHarmonicSum", -scaledHarmonic_udf("harmonic_sum", "maximum"))
    )
    mousePhenoScore = scoreAggregation.select(
        "targetFromSourceId", "scaledHarmonicSum"
    ).join(queryset, queryset.targetid == F.col("targetFromSourceId"), "right")

    return mousePhenoScore
