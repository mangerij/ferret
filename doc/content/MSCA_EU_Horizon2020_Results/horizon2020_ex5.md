!config navigation breadcrumbs=False

# H2020-MSCA-IF-2019, Example: spin wave transport

This page details how to obtain spin wave transport across the multiferroic domain boundary. In order to perform this calculation, one needs the information from two other examples (first relax the FE DW, then get the relaxed magnetic DW in the presence of the FE domain boundary). These are handled in Example 2 and 3 of the MSCA tutorials respectively.

The initial condition is a fixed $\{\mathbf{P},\mathbf{A}\} configuration with a corresponding relaxed $\{\mathbf{L},\mathbf{m}\}$ ground state. We load in the `Exodus` file from the previous example in the `Mesh` block,

This project [SCALES - 897614](https://cordis.europa.eu/project/id/897614) was funded for 2021-2023 at the [Luxembourg Institute of Science and Technology](https://www.list.lu/) under principle investigator [Jorge Íñiguez](https://sites.google.com/site/jorgeiniguezresearch/). The research was carried out within the framework of the [Marie Skłodowska-Curie Action (H2020-MSCA-IF-2019)](https://ec.europa.eu/info/funding-tenders/opportunities/portal/screen/opportunities/topic-details/msca-if-2020) fellowship.

!media media/euflag.png style=display:block;margin-left:auto;margin-right:auto;width:12%;

!content pagination previous=MSCA_EU_Horizon2020_Results/horizon2020_ex4.md
