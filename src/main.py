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

### Define queryset, what at the moment is 
### a random sample of target dataset
queryset=target.select('id').withColumnRenamed('id','targetid').sample(False, 0.3,seed=10).limit(10000)

### Define target dataset
target_path = "/Users/juanr/Desktop/Target_Engine/data_download/Parquet/target/targets/"
target = spark.read.parquet(target_path)

biotype = biotype_query(target, queryset) #### [id,biotype]
mblocation = target_location(target,biotype,collect_termsl) 
drug = drug_query(target, mblocation)

### output.write.csv()

## From questions to code
## Common: Load target and query datasets (datasets with name of targets).

### target_path = "/Users/juanr/Desktop/Target_Engine/data_download/Parquet/target/targets/"
### queryset_path = "/Users/juanr/Desktop/Target_Engine/queryset.csv"
### target = spark.read.parquet(target_path)
### queryset = spark.read.option("header", True).csv(queryset_path)


