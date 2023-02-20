from target_engine_repo.src.data_flow.target_properties_wb import (
    biotype_query,
    chemical_probes,
    clin_trials,
    constraint,
    driver_genes,
    ligand_pocket_query,
    mousemod_class,
    orthologs_mouse,
    paralogs,
    safety_query,
    target_membrane,
    tep_query,
    tissue_specific,
)
from re import X
from pyspark.sql import DataFrame, SparkSession
import pyspark.sql.functions as F
from pyspark.sql.types import StructType, StructField, StringType
import pyspark.sql.types as t

spark = (
    SparkSession.builder.master("local[*]")
    .config("spark.driver.memory", "15g")
    .appName("spark")
    .getOrCreate()
)
### Define required datasets
target_path = "/Users/juanr/Desktop/Target_Engine/data_download/Parquet/target/targets/"
target = spark.read.parquet(target_path)
interact_path = "/Users/juanr/Desktop/Target_Engine/data_download/Parquet/interaction/"
interact_db = spark.read.parquet(interact_path)
molecule_path = (
    "/Users/juanr/Desktop/Target_Engine/data_download/Parquet/drug/molecule/"
)
molecule = spark.read.parquet(molecule_path)
mouse_path = "/Users/juanr/Desktop/Target_Engine/data_download/Parquet/mousePhenotypes"
mouse = spark.read.parquet(mouse_path)
molecule_mec_path = (
    "/Users/juanr/Desktop/Target_Engine/data_download/Parquet/mechanismOfAction/"
)
molecule_mec = spark.read.parquet(molecule_mec_path)

queryset = target.select("id").withColumnRenamed("id", "targetid")
hpa_data = "/Users/juanr/Desktop/Target_Engine/proteinatlas.json"

biotype = biotype_query(target, queryset)
location = target_membrane(target, biotype)
ligand_pocket = ligand_pocket_query(target, location)
safety = safety_query(target, ligand_pocket)
selection = constraint(target, safety)
paralog = paralogs(target, selection)
ortholog = orthologs_mouse(target, paralog)
drivers = driver_genes(target, ortholog)
tep = tep_query(target, drivers)
mice = mousemod_class(mouse, tep)
chemprob = chemical_probes(target, mice)
drug_clin = clin_trials(molecule, chemprob)
tissue_specificity = tissue_specific(hpa_data, drug_clin)

info = (
    tissue_specificity
    ### Include filtering by biotype = protein coding
    ###    .filter(F.col('biotype') == 'protein_coding')
    ### Then, remove the column isProteinCoding.
    .select(
        "targetid",
        ### "Nr_biotype",
        "Nr_mb",
        "Nr_secreted",
        "Nr_Pocket",
        "Nr_Ligand",
        "Nr_Event",
        "cal_score",
        "Nr_paralogs",
        "Nr_ortholog",
        "Nr_CDG",
        "Nr_TEP",
        "Nr_Mousemodels",
        "Nr_chprob",
        "inClinicalTrials",
        "Nr_specificity",
        "Nr_distribution",
    )
    .withColumnRenamed("Nr_mb", "isInMembrane")
    .withColumnRenamed("Nr_secreted", "isSecreted")
    .withColumnRenamed("Nr_Pocket", "hasPocket")
    .withColumnRenamed("Nr_Ligand", "hasLigand")
    .withColumnRenamed("Nr_Event", "hasSafetyEvent")
    .withColumnRenamed("cal_score", "geneticConstraint")
    .withColumnRenamed("Nr_paralogs", "paralogMaxIdentityPercentage")
    .withColumnRenamed("Nr_ortholog", "mouseOrthologMaxIdentityPercentage")
    .withColumnRenamed("Nr_CDG", "isCancerDriverGene")
    .withColumnRenamed("Nr_TEP", "hasTEP")
    .withColumnRenamed("Nr_Mousemodels", "hasMouseKO")
    .withColumnRenamed("Nr_chprob", "hasHighQualityChemicalProbes")
    .withColumnRenamed("inClinicalTrials", "maxClinicalTrialPhase")
    .withColumnRenamed("Nr_specificity", "tissueSpecificity")
    .withColumnRenamed("Nr_distribution", "tissueDistribution")
)
