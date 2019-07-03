models_in <- function() return(list(ind=c("VAR","SVAR","SMVAR", "GVAR","SVECM","SVARHL",
                                          "SVARHLX","SVARHLMA", "DFM","GGVAR"),
                                    indne=c("SMVAR", "GVAR"),
                                    pop = c("MLVAR","GGVAR"),
                                    popne= c("MLVAR"),
                                    sparse=c("SVAR","SMVAR", "GVAR","SVECM","SVARHL",
                                             "SVARHLX","SVARHLMA","GGVAR"),
                                    contemporaneous = c("GVAR","GGVAR","DFM")))


