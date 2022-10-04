
def biotype_query(target, queryset):
    target_biotype = (target
    .select("id", "biotype")
    .join(queryset, target.id == queryset.targetid, "right"
    )
    .select("targetid", "biotype")
    )
    return target_biotype #### [id,biotype]

def target_location (target, queryset):
    column= (
        target
        .select('id',
            F.col('subcellularLocations'),
            F.explode_outer('subcellularLocations')
        )
        .select('id',
            F.col('col.location')
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
