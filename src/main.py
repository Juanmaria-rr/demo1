# Do not forget to activate conda environment
# Import pyspark and create session
from re import X
from pyspark.sql import DataFrame, SparkSession
import pyspark.sql.functions as F
from pyspark.sql.types import StructType, StructField, StringType

from VScode.target_engine.target.properties import biotype_query, drug_query, target_location

spark = (
    SparkSession.builder.master("local[*]")
    .config("spark.driver.memory", "15g")
    .appName("spark")
    .getOrCreate()
)

target_path = "/Users/juanr/Desktop/Target_Engine/data_download/Parquet/target/targets/"
target = spark.read.parquet(target_path)

test = biotype_query(target.select("id"), target)
test2 = target_location(test, target)
test3 = drug_query(test2, target)

test4 = target_filter_byfile(test3, path)

output.write.csv()

## From questions to code

## Common: Load target and query datasets (datasets with name of targets).

target_path = "/Users/juanr/Desktop/Target_Engine/data_download/Parquet/target/targets/"
queryset_path = "/Users/juanr/Desktop/Target_Engine/queryset.csv"
target = spark.read.parquet(target_path)
queryset = spark.read.option("header", True).csv(queryset_path)



## Daniel solution: Concatenate the actions and only with one dataframe open



### Next questions for developing code ###

# Take from the powerpoint and continue coding this.

## Q3. Where is the target expressed across the cell? Is expressed in other tissues?

## Q4. Is there any available drug for the target or other partner? 

## cargar dataset de target, interaccion y drogas

interactors="/Users/juanr/Desktop/Target_Engine/interaction/"
molecule_path= "/Users/juanr/Desktop/Target_Engine/data_download/Parquet/drug/molecule/"
target_path="/Users/juanr/Desktop/Target_Engine/data_download/Parquet/target/targets/"

interact_db=spark.read.parquet(interactors)
molecule=spark.read.parquet(molecule_path)
target=spark.read.parquet(target_path)

## target names with 100 rows
target_names=target.select('id','approvedSymbol').limit(100)

## Intact database and Score cut-off 0.42 ## Filtramos el dataset de interactors_db a >0.42



## Q5. Does it have genetic variations? Which are not tolerated?

