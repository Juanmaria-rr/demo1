
from pyspark.sql import DataFrame, SparkSession
import pyspark.sql.functions as F

spark = SparkSession.builder \
    .master('local[*]') \
    .config("spark.driver.memory", "15g") \
    .appName('spark') \
    .getOrCreate()

### Print all True in tractability_target
def target_tractab (target,queryset): 
    column=(target
    .select(
        F.col("id").alias("target_id"),
        F.explode_outer("tractability").alias("new_struct")
    )
    .select("target_id", F.col("new_struct.*"))
    # Filter for AB modality, and approval
    .filter
        (F.col("value") == "True")
    .groupBy("target_id")
    .agg(F.collect_list(F.struct(F.col('id'), F.col('value'))).alias('drug'))
    .join(queryset, queryset.targetid == F.col('target_id'),'right')
    )
    return column


### 03.09.2022
### Add drugs to parners uniquely

interact_path="/Users/juanr/Desktop/Target_Engine/data_download/Parquet/interaction/"
interact_db=spark.read.parquet(interact_path)

def partner_drugs (molecule,interact_db,queryset): 
    tar_group=(molecule
    .select(F.col('id'), F.explode(F.col('linkedTargets.rows')))
    .groupBy('col').agg(F.collect_list('id').alias('CHEMBL'))                                             
    )
    partner_drugs=(interact_db
    .filter(interact_db.sourceDatabase =='intact').select('sourceDatabase', 'targetA','targetB','scoring')
    .filter(partners.scoring > '0.42')
    .join(queryset ,queryset.targetid ==  partners_cutoff.targetA,"right")
    .dropDuplicates(['id',"targetA","targetB"])
    .join(tar_group ,F.col('targetB') == tar_group.col,"left")
    )
    return partner_drugs 


---##1## Dataset of ChEMBL drugs with linked TherapyAreas, target and disease names (individuals)
### Next function to compare linked therapeutic areas between diseases and CHEMBL


## This was done earlier, but just to remind: 

## partners_drug=dropdf_partners.join(tar_group ,dropdf_partners.targetB == tar_group.col,"left")
interact_path="/Users/juanr/Desktop/Target_Engine/data_download/Parquet/interaction/"
interact_db=spark.read.parquet(interact_path)
partners=interact_db.filter(interact_db.sourceDatabase =='intact').select('sourceDatabase', 'targetA','targetB','scoring')
partners_cutoff=partners.filter(partners.scoring > '0.42')
df_targets_partners=df_targets.join(partners_cutoff ,df_targets.id ==  partners_cutoff.targetA,"left")
dropdf_partners = df_targets_partners.dropDuplicates(['id',"targetA","targetB"])
molecule_path= "/Users/juanr/Desktop/Target_Engine/data_download/Parquet/drug/molecule/"
molecule=spark.read.parquet(molecule_path)

tojoin=molecule.select('id','linkedDiseases')

## Explode therapeutic Areas while tracking drugs

disease_therArea=(diseases.select('id','therapeuticAreas')
.select('id',F.explode_outer('therapeuticAreas').alias('IndTherArea'))
)

drug_indisease=(molecule
.select(molecule.id, F.explode(molecule.linkedTargets.rows))
.groupBy('col').agg(F.collect_list('id').alias('CHEMBL'))                                             
.join(dropdf_partners ,dropdf_partners.targetB == F.col('col'),"right")
.select(F.col('id'),F.col('targetB'),F.explode(F.col('CHEMBL')).alias('quembol'))
.groupBy(F.col('quembol').alias('drug'))
.agg(F.collect_list('id').alias('target'),F.collect_list('targetB').alias('partner'))
.join(tojoin, tojoin.id == F.col('drug'), 'left')
.select(
     'drug',
     F.col('tartget'),
     'partner',
     F.col('linkedDiseases.rows').alias('linkedDis_name'),
     F.explode_outer('linkedDis_name').alias('inddisease_name'))
.join(disease_therArea, disease_therArea.id == F.col('inddisease_name'), 'left')
.select(
     'drug',
     'inddisease_name',
     F.col('IndTherArea').alias('IndTherAreaDrug'),
     F.explode_outer(F.col('target')).alias('target'))
)

## Explode target to obtain individual ChEMBL, target, Therapy Areas and Diseases names

---##2## Dataset of targets with disease names and therapy areas (individuals) 

df2 = spark.read.option("header",True) \
     .csv("Mayas.csv")
     .select(F.col('hgnc_protein').alias('hgnc_name'),'EFO_disease','Therapy_area')

df6=(df2.join(df_targets, df_targets.hgnc_protein == df2.hgnc_name, 'left')\
.select('id','hgnc_protein','EFO_disease','Therapy_area')
.withColumn('ind_disease_therapy',F.explode(F.split('Therapy_area',',')))
.withColumn('ind_disease_name',F.explode(F.split('EFO_disease',',')))
.select('id','ind_disease_therapy','ind_disease_name')
)

### Creo que es sustituible por el queryset. 

####### df6= dataframe con target/EFO term/Therapeutic area. #######

## Anadir al df6 la columna de droga para ver la droga que tiene cada target asociada. 

## escribir el approved symbol de los target partners
partners_drug_name=(target
.select('id','approvedSymbol')
.withColumnRenamed('approvedSymbol','Symbol_targetB').withColumnRenamed('id','targetid')
.join(partners_drug, F.col('id') == partners_drug.targetB, 'right')
)
#####

partners_name=(interact_db
.filter(interact_db.sourceDatabase =='intact')
.select('sourceDatabase', 'targetA','targetB','scoring')
.filter(F.col('scoring') > '0.42')
.join(queryset ,queryset.targetid ==  F.col('targetA'),"right")
.dropDuplicates(['id',"targetA","targetB"])
.join(tar_group ,F.col('targetB') == tar_group.col,"left")
)

df7=(partners_drug_name
.select(F.col('id').alias('target'),F.explode(F.col('CHEMBL').alias('ChEMBL')))
.join(df6, df6.id == F.col('target'), 'right')
.select(
     F.col('id'),
     F.col('col').alias('ChEMBL'),
     F.col('ind_disease_therapy'),
     F.col('ind_disease_name'))
)

---##3## join the two datasets using 3 conditions: similar therapy area, similar target id and similar ChEMBL.

joined=(df7.join(drug_inddis_indtarget, (df7.ind_disease_therapy == drug_inddis_indtarget.IndTherAreaDrug) &
(df7.id == drug_inddis_indtarget.target) & (df7.ChEMBL == drug_inddis_indtarget.drug),"inner")
.select('id','ChEMBL','ind_disease_therapy','IndTherAreaDrug')
.dropDuplicates()
)

## groupBy joined3 to catch all ChEMBLs with same therapy areas
joined_grouped=joined3.groupBy(F.col('id'))\
.agg(F.collect_list(F.col('ChEMBL')).alias('CHEMBL'))

tar_group=(molecule
.select(F.col('id'), F.explode(F.col('linkedTargets.rows')))
.groupBy('col').agg(F.collect_list('id').alias('CHEMBL'))                                             
)

interact_path="/Users/juanr/Desktop/Target_Engine/data_download/Parquet/interaction/"
interact_db=spark.read.parquet(interact_path)

### Add drugs to parners uniquely

target_partner_drugs=(interact_db
.filter(interact_db.sourceDatabase =='intact').select('sourceDatabase', 'targetA','targetB','scoring')
.filter(partners.scoring > '0.42')
.join(queryset ,queryset.id ==  partners_cutoff.targetA,"right")
.dropDuplicates(['id',"targetA","targetB"])
.join(tar_group ,F.col('targetB') == tar_group.col,"left")
)
