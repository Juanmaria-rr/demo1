

def biotype_query(target, queryset):
    target_biotype = (target
    .select("id", "biotype")
    .join(queryset, target.id == queryset.targetid, "right"
    )
    .select("targetid", "biotype")
    )
    return target_biotype #### [id,biotype]

### Temporary list of membrane termSL
##  using terms we want to find. 

terms=['Cell membrane', 
'cell membrane', 
'Extracellular vesicle membrane']
regex_pattern = "|".join(terms)

collect_termsl=(target
.select(F.explode(F.col('subcellularLocations')))
.select('col.*')
.distinct()
.filter(F.col('location').rlike(regex_pattern))
.groupBy('termSL')
.agg(F.collect_list('location'))
.select('termSL')
)
collect_termsl=list(collect_termsl.select('termSL').
    toPandas()['termSL'])
regex_pattern = "|".join(collect_termsl)


def target_location (target, queryset, regex_pattern):
    #### regex_pattern= SL-XXXXX terms
    column= (
        target
        .select('id',
            F.col('subcellularLocations'),
            F.explode_outer('subcellularLocations')
        )
        .select('id',
            F.col('col.termSL'),
            F.col('col.location')
        )
        .filter(
            F.col("col.termSL").rlike(regex_pattern))
        .select('id',
            F.col('location')
            )
        .groupBy('id')
        .agg(F.collect_list('location').alias('location'))
        .join(queryset,queryset.targetid == F.col('id'),'right')
        )
    return column ### AND JOIN TO COLUMN.

def drug_query(target, queryset):
    # Getting list of tractability annotations for every target:
    tractab = (
        target
        # Exploding tractability:
        .select(
    F.col("id").alias("target_id"),
    F.explode_outer("tractability").alias("new_struct")
    )
    .select("target_id", F.col("new_struct.*"))
    # Filter for AB modality, and approval
    .filter(
        (F.col("modality") == "AB")
        & (F.col("id") == "Approved Drug")
        & (F.col("value") == "True")
    )
    # Joining with queryset:
    .join(queryset, queryset.targetid == F.col('target_id'),'right')
    )
    # Write to output:
    return tractab ### AND JOIN COLUMN TO DATASET 