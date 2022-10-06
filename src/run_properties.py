# Do not forget to activate conda environment
# Import pyspark and create session
from re import X
from pyspark.sql import DataFrame, SparkSession
import pyspark.sql.functions as F
from pyspark.sql.types import StructType, StructField, StringType

from VScode.target_engine.target.properties import biotype_query, drug_query, target_location
from target_engine_repo.src.data_flow.target_properties import chemical_probes, mousemod_class, partner_drugs

spark = (
    SparkSession.builder.master("local[*]")
    .config("spark.driver.memory", "15g")
    .appName("spark")
    .getOrCreate()
)

### Define required datasets
target_path = "/Users/juanr/Desktop/Target_Engine/data_download/Parquet/target/targets/"
target = spark.read.parquet(target_path)
interact_path="/Users/juanr/Desktop/Target_Engine/data_download/Parquet/interaction/"
interact_db=spark.read.parquet(interact_path)
molecule_path= "/Users/juanr/Desktop/Target_Engine/data_download/Parquet/drug/molecule/"
molecule=spark.read.parquet(molecule_path)
mouse_path = "/Users/juanr/Desktop/Target_Engine/data_download/Parquet/mousePhenotypes"
mouse = spark.read.parquet(mouse_path)

### Define queryset. In this case is 
### a random sample of target dataset
queryset=target.select('id').withColumnRenamed('id','targetid').sample(False, 0.3,seed=10).limit(10000)

## Concatenated functions
biotype = biotype_query(target, queryset)
location = target_location(target,biotype) 
drug = drug_query(target, location)
drug_partners = partner_drugs (molecule,interact_db,drug)
chemi_probes= chemical_probes (target,drug_partners)
mouse_models= mousemod_class (mouse,chemi_probes)

#Selection of relevant columns
info=(mouse_models
    .select(
        'targetid',
        'biotype',
        'location',
        'Approved_drugType',
        'N_partner_drug',
        'ChemicalProbes_HC',
        'Nr_mouse_models',
        'Different_PhenoClasses')
) 



