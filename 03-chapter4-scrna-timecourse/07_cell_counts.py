# Cell-count line plot before versus after QC across days and conditions.
import pandas as pd
from plotnine import (
    ggplot, aes, geom_line, geom_point, geom_text,
    labs, theme_minimal, theme, scale_color_manual,
    element_text, position_nudge, facet_wrap
)

df = pd.DataFrame({
    "day":       [1,1,1,1, 2,2,2,2, 3,3,3,3, 4,4,4,4, 5,5,5,5],
    "condition": ["MDI","MDIX"]*10,
    "qc":        ["before","before","after","after"]*5,
    "n_cells":   [4578, 4490, 4566, 4455,
                  3478, 3188, 3371, 3128,
                  2793, 2942, 2649, 2877,
                  2628, 2540, 2587, 2504,
                  574,  1629, 572,  1606],
})

df["qc"] = pd.Categorical(df["qc"], categories=["before", "after"], ordered=True)

p = (
    ggplot(df, aes(x="day", y="n_cells", color="condition", group="condition"))
    + geom_line(size=1.2)
    + geom_point(size=4)
    + geom_text(aes(label="n_cells"),
                position=position_nudge(y=200),
                size=8, show_legend=False)
    + facet_wrap("~qc")
    + scale_color_manual(values={"MDI": "#E64B35", "MDIX": "#4DBBD5"})
    + labs(
        x="Day",
        y="Number of cells",
        title="Cell number across timepoints — before vs after QC",
        color="Condition",
    )
    + theme_minimal()
    + theme(
        plot_title=element_text(size=13, weight="bold"),
        axis_title=element_text(size=11),
        strip_text=element_text(size=11, weight="bold"),
    )
)

p.save("cell_counts_qc_facet.png", width=10, height=4.5, dpi=150)
p
