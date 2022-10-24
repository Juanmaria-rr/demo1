import pyspark.sql.functions as F

def nest_prioritisations(df: DataFrame, targetcol: str) -> DataFrame:
##Nest prioritisation columns into a single column.

##Args:
##df (DataFrame): Input dataframe with prioritisation columns.
##targetcol (str): Column Name containing target ids.

##Returns:
##DataFrame: Containing target ids and `prioritisation` column

    return df.select(
    targetcol,
    F.array_except(
        F.array(
            *[
                F.when(
                    F.col(column).isNotNull(),
                    F.struct(
                        F.lit(column).alias("id"), F.col(column).alias("value")
                    ),
                )
                for column in df.columns
                if column != targetcol
            ]
        ),
        F.array(F.lit(None)),
    ).alias("prioritisations"),
)
