import pandas as pd

data = pd.read_csv("./prognostic_modules.csv", index_col=0)
cancers = pd.read_csv("./cancers.csv", index_col=0)
data.head()

key = "metabolites"

data = (
    pd.merge(data, cancers, left_on="cancer_id", right_on=cancers.index)
    .set_index("cancer_id")
    .rename(columns={"name": "cancer"})
)
group = data.groupby("cancer", as_index=False).agg(lambda x: "".join(x))
group.genes = group.genes.str.split(",")
result = group.explode("genes")
result.replace("", float("NaN"), inplace=True)
result = result.dropna(subset=["genes"])
cross = (
    pd.crosstab(
        result.cancer,
        result.genes,
        margins=True,
        values=result.cancer,
        aggfunc=pd.Series.nunique,
    )
    .fillna(0)
    .astype(int)
    .T.sort_values(by="All", ascending=False)
)

grouped_result_list = result.groupby(["genes"])["cancer"].agg(list)
r = pd.merge(grouped_result_list, cross, left_on="genes", right_on="genes")
r.cancer = r.cancer.apply(lambda x: ", ".join(x))
r = r.sort_values(by="All", ascending=False)
r.insert(1, "All", r.pop("All"))
r = r.rename(columns={"All": "count"})

uniqs = result.drop_duplicates("genes", keep=False)
uniqs = uniqs.groupby("cancer")["genes"].agg(count="count", list="unique")
uniqs["list"] = uniqs["list"].apply(lambda x: ",".join(x))
uniqs = uniqs.sort_values(by="count", ascending=False)

r.to_excel(f"results/{key}_common.xlsx")
uniqs.to_excel(f"results/{key}_unique.xlsx")

print(f" {key} - Total: {result.genes.size} - Unique: {result.genes.nunique()}")

#####################################################
#####################################################

key = "genes"
fname = "prognostic_modules"
merged = (
    pd.merge(data, cancers, left_on="cancer_id", right_on=cancers.index)
    .set_index("cancer_id")
    .rename(columns={"name": "cancer"})
)
merged[key] = merged[key].str.split(",")
result = merged.explode(key)
result = result.replace("", float("NaN"))
result = result.dropna(subset=[key])
cross = (
    pd.crosstab(
        result.cancer,
        result[key],
        margins=True,
        values=result.cancer,
        aggfunc=pd.Series.nunique,
    )
    .fillna(0)
    .astype(int)
    .T.sort_values(by="All", ascending=False)
)
grouped_result_list = result.groupby([key])["cancer"].agg(list)
r = pd.merge(grouped_result_list, cross, left_on=key, right_on=key)
r.cancer = r.cancer.apply(lambda x: ", ".join(x))
r = r.sort_values(by="All", ascending=False)
r.insert(1, "All", r.pop("All"))
r = r.rename(columns={"All": "count"})

uniqs = result.drop_duplicates(key, keep=False)
uniqs = uniqs.groupby("cancer")[key].agg(count="count", list="unique")
uniqs["list"] = uniqs["list"].apply(lambda x: ",".join(x))
uniqs = uniqs.sort_values(by="count", ascending=False)

r.to_excel(f"results/{fname}_common.xlsx")
uniqs.to_excel(f"results/{fname}_unique.xlsx")

print(f" {fname} - Total: {result[key].size} - Unique: {result[key].nunique()}")

