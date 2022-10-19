### Obtention of uniprot SLterms
path="/Users/juanr/Desktop/Target_Engine/uniprot_slterms.tsv"
df = spark.read.csv(path, sep=r'\t', header=True)
### 
df_1=(df
.select(
    'Name',
    'Subcellular location ID',
    F.explode_outer(F.split('Is a',',')).alias('Is_a_exploded'),
    F.col('Is part of'))
.select(
    'Name',
    'Subcellular location ID',
    'Is_a_exploded',
    F.explode_outer(F.split('Is part of',',')).alias('Is_part_exploded'))

## Make split and select elements seperated by ";"
# and choose the first element of the two columns
.withColumn('Is_a_exploded_SL', F.split(F.col('Is_a_exploded'), ';')
.getItem(0))
.withColumn('Is_part_SL', F.split(F.col('Is_part_exploded'), ';')
.getItem(0))
.select(F.col('Name'),F.col('Subcellular location ID').alias('SubcellID'), F.col('Is_part_SL'), F.col('Is_a_exploded_SL'))
)
##### Take parental: terms whose 'Is a" 
# column is null because they are the parents
parental=(df_1.filter(df_1.Is_a_exploded_SL.isNull())
.withColumnRenamed('Is_a_exploded_SL','Is_a')
)
child=df_1.filter(df_1.Is_a_exploded_SL.isNotNull())
child_grouped=(child
.select('SubcellID','Is_a_exploded_SL').distinct()
.groupBy('Is_a_exploded_SL')
.agg(F.collect_list(F.col('SubcellID')).alias('SubcellID_child'))
)
parental_child=(parental
.join(child_grouped, parental.SubcellID == child_grouped.Is_a_exploded_SL, 'left')
.select('Name','SubcellID','SubcellID_child')
)
### Obtention of the list with the hierarchy 
### and a new column with all the terms belonging to the name. 
hierarchy=(parental_child
.select(
    'Name',
    'SubcellID',
    'SubcellID_child',F.concat_ws(',', parental_child.SubcellID,F.col('SubcellID_child')).alias("Search"))
)
