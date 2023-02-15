source("main.r")

# Activity 2
run_all_analyses(dataset_name = "soy_2020",
                 counts_path = "../data/soy_2020/soy_2020_aggregated_counts.csv",
                 treatment_key = list(u="unexposed",e="exposed"))
run_all_analyses(dataset_name = "cor_2020",
                 counts_path = "../data/cor_2020/cor_2020_aggregated_counts.csv",
                 treatment_key = list(u="unexposed",e="exposed"))
run_all_analyses(dataset_name = "cra_2020",
                 counts_path = "../data/cra_2020/cra_2020_aggregated_counts.csv",
                 treatment_key = list(u="unexposed",e="exposed"))
run_all_analyses(dataset_name = "cas_2020",
                 counts_path = "../data/cas_2020/cas_2020_aggregated_counts.csv",
                 treatment_key = list(u="unexposed",e="exposed"))
run_all_analyses(dataset_name = "cac_2020",
                 counts_path = "../data/cac_2020/cac_2020_aggregated_counts.csv",
                 treatment_key = list(u="unexposed",e="exposed"))
run_all_analyses(dataset_name = "lbb_2021",
                 counts_path = "../data/lbb_2021/lbb_2021_aggregated_counts.csv",
                 treatment_key = list(u="unexposed",e="exposed"))
run_all_analyses(dataset_name = "cra_2021",
                 counts_path = "../data/cra_2021/cra_2021_aggregated_counts.csv",
                 treatment_key = list(u="unexposed",e="exposed"))
run_all_analyses(dataset_name = "cas_2021",
                 counts_path = "../data/cas_2021/cas_2021_aggregated_counts.csv",
                 treatment_key = list(u="unexposed",e="exposed"))
run_all_analyses(dataset_name = "cac_2021",
                 counts_path = "../data/cac_2021/cac_2021_aggregated_counts.csv",
                 treatment_key = list(u="unexposed",e="exposed"))
run_all_analyses(dataset_name = "app_2021",
                 counts_path = "../data/app_2021/app_2021_aggregated_counts.csv",
                 treatment_key = list(u="unexposed",e="exposed"))

run_all_analyses(dataset_name = "hbb_2021_t2",
                 counts_path = "../data/hbb_2021/hbb_2021_aggregated_counts.csv",
                 treatment_key = list(u="unexposed",e="exposed"),
                 filter_string = "t1|t3|t4")
run_all_analyses(dataset_name = "hbb_2020_t2",
                 counts_path = "../data/hbb_2020/hbb_2020_aggregated_counts.csv",
                 treatment_key = list(u="unexposed",e="exposed"),
                 filter_string = "t1|t3|t4")

run_all_analyses(dataset_name = "hbb_2021_t1",
                 counts_path = "../data/hbb_2021/hbb_2021_aggregated_counts.csv",
                 treatment_key = list(u="unexposed",e="exposed"),
                 filter_string = "t2|t3|t4")
run_all_analyses(dataset_name = "hbb_2021_t3",
                 counts_path = "../data/hbb_2021/hbb_2021_aggregated_counts.csv",
                 treatment_key = list(u="unexposed",e="exposed"),
                 filter_string = "t1|t2|t4")
run_all_analyses(dataset_name = "hbb_2021_t4",
                 counts_path = "../data/hbb_2021/hbb_2021_aggregated_counts.csv",
                 treatment_key = list(u="unexposed",e="exposed"),
                 filter_string = "t1|t3|t2")

run_all_analyses(dataset_name = "hbb_2020_t1",
                 counts_path = "../data/hbb_2020/hbb_2020_aggregated_counts.csv",
                 treatment_key = list(u="unexposed",e="exposed"),
                 filter_string = "t2|t3|t4")
run_all_analyses(dataset_name = "hbb_2020_t3",
                 counts_path = "../data/hbb_2020/hbb_2020_aggregated_counts.csv",
                 treatment_key = list(u="unexposed",e="exposed"),
                 filter_string = "t1|t2|t4")
run_all_analyses(dataset_name = "hbb_2020_t4",
                 counts_path = "../data/hbb_2020/hbb_2020_aggregated_counts.csv",
                 treatment_key = list(u="unexposed",e="exposed"),
                 filter_string = "t1|t3|t2")


# Activity 1
act1_meta <- read_csv("activity-1_metadata.csv")
# A_BOS
run_all_analyses(dataset_name = "a_bos_d1_low_exposure",
                 counts_path = "../data/a_bos_2021/a_bos_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d1="low_exposure"),
                 filter_string = "_d2",
                 act1_meta = act1_meta)
run_all_analyses(dataset_name = "a_bos_d2_high_exposure",
                 counts_path = "../data/a_bos_2021/a_bos_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d2="high_exposure"),
                 filter_string = "_d1",
                 act1_meta = act1_meta)
# A_FLX
run_all_analyses(dataset_name = "a_flx_d1_low_exposure",
                 counts_path = "../data/a_flx_2021/a_flx_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d1="low_exposure"),
                 filter_string = "_d2",
                 act1_meta = act1_meta)
run_all_analyses(dataset_name = "a_flx_d2_high_exposure",
                 counts_path = "../data/a_flx_2021/a_flx_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d2="high_exposure"),
                 filter_string = "_d1",
                 act1_meta = act1_meta)
# A_PYM
run_all_analyses(dataset_name = "a_pym_d1_low_exposure",
                 counts_path = "../data/a_pym_2021/a_pym_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d1="low_exposure"),
                 filter_string = "_d2",
                 act1_meta = act1_meta)
run_all_analyses(dataset_name = "a_pym_d2_high_exposure",
                 counts_path = "../data/a_pym_2021/a_pym_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d2="high_exposure"),
                 filter_string = "_d1",
                 act1_meta = act1_meta)


# B_FLY
run_all_analyses(dataset_name = "b_fly_d1_low_exposure",
                 counts_path = "../data/b_fly_2021/b_fly_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d1="low_exposure"),
                 filter_string = "_d2",
                 act1_meta = act1_meta)
run_all_analyses(dataset_name = "b_fly_d2_high_exposure",
                 counts_path = "../data/b_fly_2021/b_fly_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d2="high_exposure"),
                 filter_string = "_d1",
                 act1_meta = act1_meta)
# B_PYC
run_all_analyses(dataset_name = "b_pyc_d1_low_exposure",
                 counts_path = "../data/b_pyc_2021/b_pyc_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d1="low_exposure"),
                 filter_string = "_d2",
                 act1_meta = act1_meta)
run_all_analyses(dataset_name = "b_pyc_d2_high_exposure",
                 counts_path = "../data/b_pyc_2021/b_pyc_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d2="high_exposure"),
                 filter_string = "_d1",
                 act1_meta = act1_meta)

# c_chl
run_all_analyses(dataset_name = "c_chl_d1_low_exposure",
                 counts_path = "../data/c_chl_2021/c_chl_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d1="low_exposure"),
                 filter_string = "_d2",
                 act1_meta = act1_meta)
run_all_analyses(dataset_name = "c_chl_d2_high_exposure",
                 counts_path = "../data/c_chl_2021/c_chl_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d2="high_exposure"),
                 filter_string = "_d1",
                 act1_meta = act1_meta)
# c_spn
run_all_analyses(dataset_name = "c_spn_d1_low_exposure",
                 counts_path = "../data/c_spn_2021/c_spn_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d1="low_exposure"),
                 filter_string = "_d2",
                 act1_meta = act1_meta)
run_all_analyses(dataset_name = "c_spn_d2_high_exposure",
                 counts_path = "../data/c_spn_2021/c_spn_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d2="high_exposure"),
                 filter_string = "_d1",
                 act1_meta = act1_meta)
# c_spr
run_all_analyses(dataset_name = "c_spr_d1_low_exposure",
                 counts_path = "../data/c_spr_2021/c_spr_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d1="low_exposure"),
                 filter_string = "_d2",
                 act1_meta = act1_meta)
run_all_analyses(dataset_name = "c_spr_d2_high_exposure",
                 counts_path = "../data/c_spr_2021/c_spr_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d2="high_exposure"),
                 filter_string = "_d1",
                 act1_meta = act1_meta)

# d_flp
run_all_analyses(dataset_name = "d_flp_d1_low_exposure",
                 counts_path = "../data/d_flp_2021/d_flp_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d1="low_exposure"),
                 filter_string = "_d2",
                 act1_meta = act1_meta)
run_all_analyses(dataset_name = "d_flp_d2_high_exposure",
                 counts_path = "../data/d_flp_2021/d_flp_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d2="high_exposure"),
                 filter_string = "_d1",
                 act1_meta = act1_meta)
# d_sul
run_all_analyses(dataset_name = "d_sul_d1_low_exposure",
                 counts_path = "../data/d_sul_2021/d_sul_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d1="low_exposure"),
                 filter_string = "_d2",
                 act1_meta = act1_meta)
run_all_analyses(dataset_name = "d_sul_d2_high_exposure",
                 counts_path = "../data/d_sul_2021/d_sul_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d2="high_exposure"),
                 filter_string = "_d1",
                 act1_meta = act1_meta)

# e_gly
run_all_analyses(dataset_name = "e_gly_d1_low_exposure",
                 counts_path = "../data/e_gly_2021/e_gly_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d1="low_exposure"),
                 filter_string = "_d2",
                 act1_meta = act1_meta)
run_all_analyses(dataset_name = "e_gly_d2_high_exposure",
                 counts_path = "../data/e_gly_2021/e_gly_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d2="high_exposure"),
                 filter_string = "_d1",
                 act1_meta = act1_meta)
# e_met
run_all_analyses(dataset_name = "e_met_d1_low_exposure",
                 counts_path = "../data/e_met_2021/e_met_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d1="low_exposure"),
                 filter_string = "_d2",
                 act1_meta = act1_meta)
run_all_analyses(dataset_name = "e_met_d2_high_exposure",
                 counts_path = "../data/e_met_2021/e_met_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d2="high_exposure"),
                 filter_string = "_d1",
                 act1_meta = act1_meta)

# CLO
run_all_analyses(dataset_name = "clo_2020_d1_sublethal",
                 counts_path = "../data/clo_2020/clo_2020_aggregated_counts.csv",
                 treatment_key = list(d0="control",d1="sublethal"),
                 filter_string = "_d2",
                 act1_meta = act1_meta)
run_all_analyses(dataset_name = "clo_2020_d2_acute",
                 counts_path = "../data/clo_2020/clo_2020_aggregated_counts.csv",
                 treatment_key = list(d0="control",d2="acute"),
                 filter_string = "_d1",
                 act1_meta = act1_meta)
# THI
run_all_analyses(dataset_name = "thi_2020_d1_sublethal",
                 counts_path = "../data/thi_2020/thi_2020_aggregated_counts.csv",
                 treatment_key = list(d0="control",d1="sublethal"),
                 filter_string = "_d2",
                 act1_meta = act1_meta)
run_all_analyses(dataset_name = "thi_2020_d2_acute",
                 counts_path = "../data/thi_2020/thi_2020_aggregated_counts.csv",
                 treatment_key = list(d0="control",d2="acute"),
                 filter_string = "_d1",
                 act1_meta = act1_meta)
# CTX
run_all_analyses(dataset_name = "ctx_2020_dC_clothianidin",
                 counts_path = "../data/ctx_2020/ctx_dC_2020_aggregated_counts.csv",
                 treatment_key = list(d0="control",dC="sublethal"),
                 act1_meta = act1_meta)
run_all_analyses(dataset_name = "ctx_2020_dT_thiamethoxam",
                 counts_path = "../data/ctx_2020/ctx_dT_2020_aggregated_counts.csv",
                 treatment_key = list(d0="control",dT="acute"),
                 act1_meta = act1_meta)
# CFS
run_all_analyses(dataset_name = "cfs_2021_dC_chlorantraniliprole",
                 counts_path = "../data/cfs_dC_2021/cfs_dC_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",dC="chlorantraniliprole"),
                 act1_meta = act1_meta)
run_all_analyses(dataset_name = "cfs_2021_dF_flupyradifurone",
                 counts_path = "../data/cfs_dF_2021/cfs_dF_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",dF="flupyradifurone"),
                 act1_meta = act1_meta)
run_all_analyses(dataset_name = "cfs_2021_dS_sulfoxaflor",
                 counts_path = "../data/cfs_dS_2021/cfs_dS_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",dS="sulfoxaflor"),
                 act1_meta = act1_meta)

# AFB
run_all_analyses(dataset_name = "afb_2021_d1_subclinical",
                 counts_path = "../data/afb_2021/afb_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d1="subclinical"),
                 filter_string = "_d2",
                 act1_meta = act1_meta)
run_all_analyses(dataset_name = "afb_2021_d2_clinical",
                 counts_path = "../data/afb_2021/afb_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d2="clinical"),
                 filter_string = "_d1",
                 act1_meta = act1_meta)

# IAP
run_all_analyses(dataset_name = "iap_2021",
                 counts_path = "../data/iap_2021/iap_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d1="infected"),
                 act1_meta = act1_meta)

# NOS
run_all_analyses(dataset_name = "nos_2021_d1_low_dose",
                 counts_path = "../data/nos_2021/nos_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d1="low_dose"),
                 filter_string = "_d2",
                 act1_meta = act1_meta)
run_all_analyses(dataset_name = "nos_2021_d2_high_dose",
                 counts_path = "../data/nos_2021/nos_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d2="high_dose"),
                 filter_string = "_d1",
                 act1_meta = act1_meta)

# VAR
run_all_analyses(dataset_name = "var_2021_d1_low_dose",
                 counts_path = "../data/var_2021/var_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d1="low_dose"),
                 filter_string = "_d2",
                 act1_meta = act1_meta)
run_all_analyses(dataset_name = "var_2021_d2_high_dose",
                 counts_path = "../data/var_2021/var_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d2="high_dose"),
                 filter_string = "_d1",
                 act1_meta = act1_meta)
# OXY
run_all_analyses(dataset_name = "oxy_2021",
                 counts_path = "../data/oxy_2021/oxy_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d1="oxytetracycline"),
                 act1_meta = act1_meta)

# PDV
run_all_analyses(dataset_name = "pdv_2021_d1_nutritious_monofloral",
                 counts_path = "../data/pdv_2021/pdv_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d1="nutritious_monofloral"),
                 filter_string = "_d2|_t1",
                 act1_meta = act1_meta)
run_all_analyses(dataset_name = "pdv_2021_d2_non_nutritious_polyfloral",
                 counts_path = "../data/pdv_2021/pdv_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d2="non_nutritious_polyfloral"),
                 filter_string = "_d1|_t1",
                 act1_meta = act1_meta)

# PRE
run_all_analyses(dataset_name = "pre_2021_d1_restricted",
                 counts_path = "../data/pre_2021/pre_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d1="restricted"),
                 filter_string = "_d2|_t1",
                 act1_meta = act1_meta)
run_all_analyses(dataset_name = "pre_2021_d2_supplemented",
                 counts_path = "../data/pre_2021/pre_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d2="supplemented"),
                 filter_string = "_d1|_t1",
                 act1_meta = act1_meta)

# AMZ
run_all_analyses(dataset_name = "amz_2021",
                 counts_path = "../data/amz_2021/amz_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d1="amitraz"),
                 act1_meta = act1_meta)
# OXA
run_all_analyses(dataset_name = "oxa_2021",
                 counts_path = "../data/oxa_2021/oxa_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d1="oxalic_acid"),
                 act1_meta = act1_meta)
# CHA
run_all_analyses(dataset_name = "cha_2021",
                 counts_path = "../data/cha_2021/cha_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d1="chalkbrood"),
                 act1_meta = act1_meta)

# FMA
run_all_analyses(dataset_name = "fma_2021_d1_low_dose",
                 counts_path = "../data/fma_2021/fma_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d1="low_dose"),
                 filter_string = "_d2",
                 act1_meta = act1_meta)
run_all_analyses(dataset_name = "fma_2021_d2_high_dose",
                 counts_path = "../data/fma_2021/fma_2021_aggregated_counts.csv",
                 treatment_key = list(d0="control",d2="high_dose"),
                 filter_string = "_d1",
                 act1_meta = act1_meta)
