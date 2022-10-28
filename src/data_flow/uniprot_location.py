## modification of parental dataset (we do not filter by null)
### This seem to be the reliable part to work with the SLterms. 
path="/Users/juanr/Desktop/Target_Engine/uniprot_slterms.tsv"
df = spark.read.csv(path, sep=r'\t', header=True)
### 
## Process the data to obtain the correct values of the separeted by ;
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
######

parental=(df_1
.select('Name','SubcellID','Is_a_exploded_SL').distinct()
.withColumnRenamed('Is_a_exploded_SL','Is_a')
)
## modification of child (we do not filter by not null)
child=(df_1
.select('SubcellID','Is_a_exploded_SL').distinct()
.groupBy('Is_a_exploded_SL')
.agg(F.collect_list(F.col('SubcellID')).alias('SubcellID_child'))
)

parental_child=(parental
.join(child, parental.SubcellID == child.Is_a_exploded_SL, 'left')
.select('Name','SubcellID','Is_a','SubcellID_child')
)

cousins=(df_1
.groupBy(F.col('Is_part_SL'))
.agg(F.collect_list('SubcellID').alias('SubcellID_are_part'))
)
## we join with cousins: 
### unimos los cousins
parent_child_cousins=(parental_child
.join(cousins, cousins.Is_part_SL == parental_child.SubcellID, 'left')
.select(
    F.col('Name'),
    F.col('SubcellID'),
    F.col('Is_a'),
    F.col('SubcellID_child').alias('Child_SLterms/parent_of'),
    F.col('SubcellID_are_part').alias('Contains_SLterms'),
)
.withColumn('concat', 
    F.concat_ws(',',
        F.col('SubcellID'),
        F.col('Child_SLterms/parent_of'),
        F.col('Contains_SLterms'))
    )
.withColumn('toSearch', F.split(
    F.col('concat'),",")
    )
)
membrane_terms=(parent_child_cousins
.filter(F.col('Name') == 'Cell membrane')
.select(parent_child_cousins.toSearch).rdd.flatMap(lambda x: x).collect()[0]
)
