''' 
workspace for creating mech_dicts for use in milstein neurons
'''
import pickle

ca1_milstein_1 = {
            # soma
            #------------------------------
            'soma': {
                'cable':{
                    'Ra':{
                        'value':150.
                    },
                    'cm':{
                        'value':1.
                    },
                },
                'ions':{
                    'ek':{
                        'value':-77.
                    },
                },
                'h_milstein':{
                    'eh':{
                        'value':-30.
                    },
                    'ghbar':{
                        'value':2.43e-9
                    },
                    'vhalfl':{
                        'value':-73.
                    },
                },
                'kap':{
                    'gkabar':{
                        'value':.0305
                    },
                },
                'kdr':{
                    'gkdrbar':{
                        'value':.0478
                    },
                },
                'km2':{
                    'gkmbar':{
                        'value':.0015
                    },
                },
                'nas':{
                    'ar':{
                        'value':1.
                    },
                    'gbar':{
                        'value':.04
                    },
                    'sh':{
                        'value':1.7
                    },
                },
                'pas':{
                    'e':{
                        'value':-63.
                    },
                    'g':{
                        'value':1.52e-6
                    },
                },
            },
            # axon
            #-----------------------------------
            'axon': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'ions':{
                    'ek':{
                        'origin':'soma'
                    },
                },
                'kap':{
                    'gkabar':{
                        'value':.090585
                    },
                },
                'kdr':{
                    'gkdrbar':{
                        'origin':'soma'
                    },
                },
                'km2':{
                    'gkmbar':{
                        'origin':'ais'
                    },
                },
                'nax_milstein':{
                    'gbar':{
                        'value':.0936
                    },
                    'sh':{
                        'origin':'axon_hill'
                    },
                },
                'pas':{
                    'e':{
                        'value':-67.
                    },
                    'g':{
                        'origin':'soma'
                    },
                },
            },
            # ais
            #-------------------------------------
            'ais': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'ions':{
                    'ek':{
                        'origin':'soma'
                    },
                },
                'kap':{
                    'gkabar':{
                        'value':.090585
                    },
                },
                'kdr':{
                    'gkdrbar':{
                        'origin':'soma'
                    },
                },
                'km2':{
                    'gkmbar':{
                        'value':.0075
                    },
                },
                'nax_milstein':{
                    'gbar':{
                        'value':.41184
                    },
                    'sh':{
                        'origin':'axon_hill'
                    },
                    'sha':{
                        'value':-3.6
                    },
                },
                'pas':{
                    'e':{
                        'value':-67.
                    },
                    'g':{
                        'origin':'soma'
                    },
                },
            },
            # axon_hill
            #---------------------------------------
            'axon_hill': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'ions':{
                    'ek':{
                        'origin':'soma'
                    },
                },
                'kap':{
                    'gkabar':{
                        'origin':'soma'
                    },
                },
                'kdr':{
                    'gkdrbar':{
                        'origin':'soma'
                    },
                },
                'km2':{
                    'gkmbar':{
                        'origin':'soma'
                    },
                },
                'nax_milstein':{
                    'gbar':{
                        'value':.04
                    },
                    'sh':{
                        'value':1.7
                    },
                },
                'pas':{
                    'e':{
                        'value':-67.
                    },
                    'g':{
                        'origin':'soma'
                    },
                },
            },
            # basal
            #-----------------------------------------
            'basal': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'ions':{
                    'ek':{
                        'origin':'soma'
                    },
                },
                'h_milstein':{
                    'eh':{
                        'origin':'soma'
                    },
                    'ghbar':{
                        'origin':'soma'
                    },
                    'vhalfl':[
                        {'max_loc':75., 'origin':'soma'},
                        {'min_loc':75., 'origin':'soma', 'value':-81}
                    ],
                },
                'kad':{
                    # 'eh':{
                    #     'origin':'soma'
                    # },
                    'gkabar':[
                        {'max_loc':75., 'origin':'soma', 'value':0.},
                        {'min_loc':75., 'max_loc':300., 'origin':'soma', 'slope':.000183, 'value':0.044225},
                        {'min_loc':300., 'origin':'soma', 'value':0.0854},
                    ],
                },
                'kap':{
                    'gkabar':[
                        {'min_loc':75., 'origin':'soma', 'value':0.},
                        {'max_loc':75., 'origin':'soma', 'slope':0.000183},
                    ],
                },
                'kdr':{
                    'gkdrbar':{
                        'origin':'soma'
                    },
                },
                'nas':{
                    'ar':{
                        'value':0.8
                    },
                    'gbar':{
                        'min':0., 'origin':'soma', 'slope':-0.0002,
                    },
                    'sh':{
                        'value':5.
                    },
                },
                'pas':{
                    'e':{
                        'value':-67.
                    },
                    'g':{
                        'origin':'soma'
                    },
                },
                'synapse':{
                    'AMPA_KIN':{
                            'gmax':{
                                'value':0.00079
                            }
                    },
                    'NMDA_KIN5':{
                            'gmax':{
                                'value':0.003026
                            }
                    },
                },
            },
            # trunk
            #------------------------------------------------
            'trunk': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'ions':{
                    'ek':{
                        'origin':'soma'
                    },
                },
                'h_milstein':{
                    'eh':{
                        'origin':'soma'
                    },
                    'ghbar':{
                        'origin':'soma',
                        'slope':0.00185,
                        'tau':383.1,
                        'xhalf':260.7,
                    },
                    'vhalfl':[
                        {'max_loc':75., 'origin':'soma'},
                        {'min_loc':75., 'origin':'soma', 'value':-81}
                    ],
                },
                'kad':{
                    # 'eh':{
                    #     'origin':'soma'
                    # },
                    'gkabar':[
                        {'max_loc':75., 'origin':'soma', 'value':0.},
                        {'min_loc':75., 'max_loc':300., 'origin':'soma', 'slope':.000183, 'value':0.044225},
                        {'min_loc':300., 'origin':'soma', 'value':0.0854},
                    ],
                },
                'kap':{
                    'gkabar':[
                        {'min_loc':75., 'origin':'soma', 'value':0.},
                        {'max_loc':75., 'origin':'soma', 'slope':0.000183},
                    ],
                },
                'kdr':{
                    'gkdrbar':{
                        'origin':'soma'
                    },
                },
                'nas':{
                    'ar':{
                        'value':0.8
                    },
                    'gbar':{
                        'origin':'soma',
                    },
                    'sh':{
                        'origin':'soma',
                    },
                    'sha':{
                        'value':5.
                    },
                },
                'pas':{
                    'e':{
                        'origin':'soma',
                        'slope':-1.582714325,
                        'tau':121.9,
                    },
                    'g':{
                        'origin':'soma',
                        'slope':-1.63e-5,
                        'tau':121.9,
                    },
                },
                'synapse':{
                    'AMPA_KIN':{
                            'gmax':{
                                'origin':'soma',
                                'slope':0.000297,
                                'tau':132.8,
                                'value':0.00079,
                            },
                    },
                    'NMDA_KIN5':{
                            'gmax':{
                                'value':0.003026
                            }
                    },
                },
            },
            # apical
            #--------------------------------------------------------
            'apical': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'ions':{
                    'ek':{
                        'origin':'soma'
                    },
                },
                'h_milstein':{
                    'eh':{
                        'origin':'soma'
                    },
                    'ghbar':{
                        'origin':'trunk',
                    },
                    'vhalfl':[
                        {'max_loc':75., 'origin':'soma'},
                        {'min_loc':75., 'origin':'soma', 'value':-81}
                    ],
                },
                'kad':{
                    # 'eh':{
                    #     'origin':'soma'
                    # },
                    'gkabar':[
                        {'max_loc':75., 'origin':'soma', 'value':0.},
                        {'min_loc':75., 'max_loc':300., 'origin':'soma', 'slope':.000183, 'value':0.044225},
                        {'min_loc':300., 'origin':'soma', 'value':0.0854},
                    ],
                },
                'kap':{
                    'gkabar':[
                        {'min_loc':75., 'origin':'soma', 'value':0.},
                        {'max_loc':75., 'origin':'soma', 'slope':0.000183},
                    ],
                },
                'kdr':{
                    'gkdrbar':{
                        'origin':'soma'
                    },
                },
                'nas':{
                    'ar':{
                        'value':0.8
                    },
                    'gbar':{
                        'origin':'trunk',
                        'min':0.0,
                        'slope':-0.0002,
                    },
                    'sh':{
                        'origin':'soma',
                    },
                    'sha':{
                        'value':5.
                    },
                },
                'pas':{
                    'e':{
                        'origin':'trunk',
                    },
                    'g':{
                        'origin':'trunk',
                    },
                },
                'synapse':{
                    'AMPA_KIN':{
                            'gmax':{
                                'origin':'trunk',
                            },
                    },
                    'NMDA_KIN5':{
                            'gmax':{
                                'origin':'trunk',
                            }
                    },
                },
            },
            # tuft
            #---------------------------------------------
            'tuft': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'ions':{
                    'ek':{
                        'origin':'soma'
                    },
                },
                'h_milstein':{
                    'eh':{
                        'origin':'soma'
                    },
                    'ghbar':{
                        'origin':'trunk',
                    },
                    'vhalfl':[
                        {'max_loc':75., 'origin':'soma'},
                        {'min_loc':75., 'origin':'soma', 'value':-81}
                    ],
                },
                'kad':{
                    # 'eh':{
                    #     'origin':'soma'
                    # },
                    'gkabar':[
                        {'max_loc':75., 'origin':'soma', 'value':0.},
                        {'min_loc':75., 'max_loc':300., 'origin':'soma', 'slope':.000183, 'value':0.044225},
                        {'min_loc':300., 'origin':'soma', 'value':0.0854},
                    ],
                },
                'kap':{
                    'gkabar':[
                        {'min_loc':75., 'origin':'soma', 'value':0.},
                        {'max_loc':75., 'origin':'soma', 'slope':0.000183},
                    ],
                },
                'kdr':{
                    'gkdrbar':{
                        'origin':'soma'
                    },
                },
                'nas':{
                    'ar':{
                        'value':0.8
                    },
                    'gbar':{
                        'origin':'trunk',
                    },
                    'sh':{
                        'origin':'soma',
                    },
                    'sha':{
                        'value':5.
                    },
                },
                'pas':{
                    'e':{
                        'origin':'trunk',
                    },
                    'g':{
                        'origin':'trunk',
                    },
                },
                'synapse':{
                    'AMPA_KIN':{
                            'gmax':{
                                'origin':'trunk',
                            },
                    },
                    'NMDA_KIN5':{
                            'gmax':{
                                'origin':'trunk',
                            }
                    },
                },
            },
            # spine_head
            #-------------------------------------------
            'spine_head': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'pas':{
                    'e':{
                        'value':-67,
                    },
                    'g':{
                        'origin':'parent',
                    },
                },
            },
            # spine_neck
            #----------------------------------------------
            'spine_neck': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'pas':{
                    'e':{
                        'value':-67,
                    },
                    'g':{
                        'origin':'parent',
                    },
                },
            },
}

ca1_milstein_high_g = {
            # soma
            #------------------------------
            'soma': {
                'cable':{
                    'Ra':{
                        'value':150.
                    },
                    'cm':{
                        'value':1.
                    },
                },
                'ions':{
                    'ek':{
                        'value':-77.
                    },
                },
                'h_milstein':{
                    'eh':{
                        'value':-30.
                    },
                    'ghbar':{
                        'value':2.43e-9
                    },
                    'vhalfl':{
                        'value':-73.
                    },
                },
                'kap':{
                    'gkabar':{
                        'value':.0305
                    },
                },
                'kdr':{
                    'gkdrbar':{
                        'value':.0478
                    },
                },
                'km2':{
                    'gkmbar':{
                        'value':.0015
                    },
                },
                'nas':{
                    'ar':{
                        'value':1.
                    },
                    'gbar':{
                        'value':0.*.04
                    },
                    'sh':{
                        'value':1.7
                    },
                },
                'pas':{
                    'e':{
                        'value':-63.
                    },
                    'g':{
                        'value':3.5e-5#1.52e-6
                    },
                },
            },
            # axon
            #-----------------------------------
            'axon': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'ions':{
                    'ek':{
                        'origin':'soma'
                    },
                },
                'kap':{
                    'gkabar':{
                        'value':.090585
                    },
                },
                'kdr':{
                    'gkdrbar':{
                        'origin':'soma'
                    },
                },
                'km2':{
                    'gkmbar':{
                        'origin':'ais'
                    },
                },
                'nax_milstein':{
                    'gbar':{
                        'value':.0936
                    },
                    'sh':{
                        'origin':'axon_hill'
                    },
                },
                'pas':{
                    'e':{
                        'value':-67.
                    },
                    'g':{
                        'origin':'soma'
                    },
                },
            },
            # ais
            #-------------------------------------
            'ais': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'ions':{
                    'ek':{
                        'origin':'soma'
                    },
                },
                'kap':{
                    'gkabar':{
                        'value':.090585
                    },
                },
                'kdr':{
                    'gkdrbar':{
                        'origin':'soma'
                    },
                },
                'km2':{
                    'gkmbar':{
                        'value':.0075
                    },
                },
                'nax_milstein':{
                    'gbar':{
                        'value':.41184
                    },
                    'sh':{
                        'origin':'axon_hill'
                    },
                    'sha':{
                        'value':-3.6
                    },
                },
                'pas':{
                    'e':{
                        'value':-67.
                    },
                    'g':{
                        'origin':'soma'
                    },
                },
            },
            # axon_hill
            #---------------------------------------
            'axon_hill': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'ions':{
                    'ek':{
                        'origin':'soma'
                    },
                },
                'kap':{
                    'gkabar':{
                        'origin':'soma'
                    },
                },
                'kdr':{
                    'gkdrbar':{
                        'origin':'soma'
                    },
                },
                'km2':{
                    'gkmbar':{
                        'origin':'soma'
                    },
                },
                'nax_milstein':{
                    'gbar':{
                        'value':.04
                    },
                    'sh':{
                        'value':1.7
                    },
                },
                'pas':{
                    'e':{
                        'value':-67.
                    },
                    'g':{
                        'origin':'soma'
                    },
                },
            },
            # basal
            #-----------------------------------------
            'basal': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'ions':{
                    'ek':{
                        'origin':'soma'
                    },
                },
                'h_milstein':{
                    'eh':{
                        'origin':'soma'
                    },
                    'ghbar':{
                        'origin':'soma'
                    },
                    'vhalfl':[
                        {'max_loc':75., 'origin':'soma'},
                        {'min_loc':75., 'origin':'soma', 'value':-81}
                    ],
                },
                'kad':{
                    # 'eh':{
                    #     'origin':'soma'
                    # },
                    'gkabar':[
                        {'max_loc':75., 'origin':'soma', 'value':0.},
                        {'min_loc':75., 'max_loc':300., 'origin':'soma', 'slope':.000183, 'value':0.044225},
                        {'min_loc':300., 'origin':'soma', 'value':0.0854},
                    ],
                },
                'kap':{
                    'gkabar':[
                        {'min_loc':75., 'origin':'soma', 'value':0.},
                        {'max_loc':75., 'origin':'soma', 'slope':0.000183},
                    ],
                },
                'kdr':{
                    'gkdrbar':{
                        'origin':'soma'
                    },
                },
                'nas':{
                    'ar':{
                        'value':0.8
                    },
                    'gbar':{
                        'min':0., 'origin':'soma', 'slope':-0.0002,
                    },
                    'sh':{
                        'value':5.
                    },
                },
                'pas':{
                    'e':{
                        'value':-67.
                    },
                    'g':{
                        'origin':'soma'
                    },
                },
                'synapse':{
                    'AMPA_KIN':{
                            'gmax':{
                                'value':0.00079
                            }
                    },
                    'NMDA_KIN5':{
                            'gmax':{
                                'value':0.003026
                            }
                    },
                },
            },
            # trunk
            #------------------------------------------------
            'trunk': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'ions':{
                    'ek':{
                        'origin':'soma'
                    },
                },
                'h_milstein':{
                    'eh':{
                        'origin':'soma'
                    },
                    'ghbar':{
                        'origin':'soma',
                        'slope':0.00185,
                        'tau':383.1,
                        'xhalf':260.7,
                    },
                    'vhalfl':[
                        {'max_loc':75., 'origin':'soma'},
                        {'min_loc':75., 'origin':'soma', 'value':-81}
                    ],
                },
                'kad':{
                    # 'eh':{
                    #     'origin':'soma'
                    # },
                    'gkabar':[
                        {'max_loc':75., 'origin':'soma', 'value':0.},
                        {'min_loc':75., 'max_loc':300., 'origin':'soma', 'slope':.000183, 'value':0.044225},
                        {'min_loc':300., 'origin':'soma', 'value':0.0854},
                    ],
                },
                'kap':{
                    'gkabar':[
                        {'min_loc':75., 'origin':'soma', 'value':0.},
                        {'max_loc':75., 'origin':'soma', 'slope':0.000183},
                    ],
                },
                'kdr':{
                    'gkdrbar':{
                        'origin':'soma'
                    },
                },
                'nas':{
                    'ar':{
                        'value':0.8
                    },
                    'gbar':{
                        'origin':'soma',
                    },
                    'sh':{
                        'origin':'soma',
                    },
                    'sha':{
                        'value':5.
                    },
                },
                'pas':{
                    'e':{
                        'origin':'soma',
                        'slope':-1.582714325,
                        'tau':121.9,
                    },
                    'g':{
                        'origin':'soma',
                        'slope':-1.63e-5,
                        'tau':121.9,
                    },
                },
                'synapse':{
                    'AMPA_KIN':{
                            'gmax':{
                                'origin':'soma',
                                'slope':0.000297,
                                'tau':132.8,
                                'value':0.00079,
                            },
                    },
                    'NMDA_KIN5':{
                            'gmax':{
                                'value':0.003026
                            }
                    },
                },
            },
            # apical
            #--------------------------------------------------------
            'apical': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'ions':{
                    'ek':{
                        'origin':'soma'
                    },
                },
                'h_milstein':{
                    'eh':{
                        'origin':'soma'
                    },
                    'ghbar':{
                        'origin':'trunk',
                    },
                    'vhalfl':[
                        {'max_loc':75., 'origin':'soma'},
                        {'min_loc':75., 'origin':'soma', 'value':-81}
                    ],
                },
                'kad':{
                    # 'eh':{
                    #     'origin':'soma'
                    # },
                    'gkabar':[
                        {'max_loc':75., 'origin':'soma', 'value':0.},
                        {'min_loc':75., 'max_loc':300., 'origin':'soma', 'slope':.000183, 'value':0.044225},
                        {'min_loc':300., 'origin':'soma', 'value':0.0854},
                    ],
                },
                'kap':{
                    'gkabar':[
                        {'min_loc':75., 'origin':'soma', 'value':0.},
                        {'max_loc':75., 'origin':'soma', 'slope':0.000183},
                    ],
                },
                'kdr':{
                    'gkdrbar':{
                        'origin':'soma'
                    },
                },
                'nas':{
                    'ar':{
                        'value':0.8
                    },
                    'gbar':{
                        'origin':'trunk',
                        'min':0.0,
                        'slope':-0.0002,
                    },
                    'sh':{
                        'origin':'soma',
                    },
                    'sha':{
                        'value':5.
                    },
                },
                'pas':{
                    'e':{
                        'origin':'trunk',
                    },
                    'g':{
                        'origin':'trunk',
                    },
                },
                'synapse':{
                    'AMPA_KIN':{
                            'gmax':{
                                'origin':'trunk',
                            },
                    },
                    'NMDA_KIN5':{
                            'gmax':{
                                'origin':'trunk',
                            }
                    },
                },
            },
            # tuft
            #---------------------------------------------
            'tuft': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'ions':{
                    'ek':{
                        'origin':'soma'
                    },
                },
                'h_milstein':{
                    'eh':{
                        'origin':'soma'
                    },
                    'ghbar':{
                        'origin':'trunk',
                    },
                    'vhalfl':[
                        {'max_loc':75., 'origin':'soma'},
                        {'min_loc':75., 'origin':'soma', 'value':-81}
                    ],
                },
                'kad':{
                    # 'eh':{
                    #     'origin':'soma'
                    # },
                    'gkabar':[
                        {'max_loc':75., 'origin':'soma', 'value':0.},
                        {'min_loc':75., 'max_loc':300., 'origin':'soma', 'slope':.000183, 'value':0.044225},
                        {'min_loc':300., 'origin':'soma', 'value':0.0854},
                    ],
                },
                'kap':{
                    'gkabar':[
                        {'min_loc':75., 'origin':'soma', 'value':0.},
                        {'max_loc':75., 'origin':'soma', 'slope':0.000183},
                    ],
                },
                'kdr':{
                    'gkdrbar':{
                        'origin':'soma'
                    },
                },
                'nas':{
                    'ar':{
                        'value':0.8
                    },
                    'gbar':{
                        'origin':'trunk',
                    },
                    'sh':{
                        'origin':'soma',
                    },
                    'sha':{
                        'value':5.
                    },
                },
                'pas':{
                    'e':{
                        'origin':'trunk',
                    },
                    'g':{
                        'origin':'trunk',
                    },
                },
                'synapse':{
                    'AMPA_KIN':{
                            'gmax':{
                                'origin':'trunk',
                            },
                    },
                    'NMDA_KIN5':{
                            'gmax':{
                                'origin':'trunk',
                            }
                    },
                },
            },
            # spine_head
            #-------------------------------------------
            'spine_head': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'pas':{
                    'e':{
                        'value':-67,
                    },
                    'g':{
                        'origin':'parent',
                    },
                },
            },
            # spine_neck
            #----------------------------------------------
            'spine_neck': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'pas':{
                    'e':{
                        'value':-67,
                    },
                    'g':{
                        'origin':'parent',
                    },
                },
            },
}
ca1_milstein_high_g_2 = {
            #soma
            #------------------------------
            'soma': {
                'cable':{
                    'Ra':{
                        'value':150.
                    },
                    'cm':{
                        'value':1.
                    },
                },
                'ions':{
                    'ek':{
                        'value':-77.
                    },
                },
                'h_milstein':{
                    'eh':{
                        'value':-30.
                    },
                    'ghbar':{
                        'value':2.43e-9
                    },
                    'vhalfl':{
                        'value':-73.
                    },
                },
                'kap':{
                    'gkabar':{
                        'value':.0305
                    },
                },
                'kdr':{
                    'gkdrbar':{
                        # 'value':.0478
                        'value':4*.0478
                    },
                },
                'km2':{
                    'gkmbar':{
                        'value':.0015
                    },
                },
                'nas':{
                    'ar':{
                        'value':1.
                    },
                    'gbar':{
                        'value':0.5*.04
                    },
                    'sh':{
                        'value':1.7
                    },
                },
                'pas':{
                    'e':{
                        'value':-67.
                    },
                    'g':{
                        # 'value':1.52e-6
                        'value':3.5e-5
                    },
                },
            },
            # axon
            #-----------------------------------
            'axon': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'ions':{
                    'ek':{
                        'origin':'soma'
                    },
                },
                'kap':{
                    'gkabar':{
                        'value':.090585
                    },
                },
                'kdr':{
                    'gkdrbar':{
                        'origin':'soma'
                    },
                },
                'km2':{
                    'gkmbar':{
                        'origin':'ais'
                    },
                },
                'nax_milstein':{
                    'gbar':{
                        'value':0.5*.0936
                    },
                    'sh':{
                        'origin':'axon_hill'
                    },
                },
                'pas':{
                    'e':{
                        'value':-67.
                    },
                    'g':{
                        'origin':'soma'
                    },
                },
            },
            # ais
            #-------------------------------------
            'ais': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'ions':{
                    'ek':{
                        'origin':'soma'
                    },
                },
                'kap':{
                    'gkabar':{
                        'value':.090585
                    },
                },
                'kdr':{
                    'gkdrbar':{
                        'origin':'soma'
                    },
                },
                'km2':{
                    'gkmbar':{
                        'value':.0075
                    },
                },
                'nax_milstein':{
                    'gbar':{
                        'value':0.5*.41184
                    },
                    'sh':{
                        'origin':'axon_hill'
                    },
                    'sha':{
                        'value':-3.6
                    },
                },
                'pas':{
                    'e':{
                        'value':-67.
                    },
                    'g':{
                        'origin':'soma'
                    },
                },
            },
            # axon_hill
            #---------------------------------------
            'axon_hill': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'ions':{
                    'ek':{
                        'origin':'soma'
                    },
                },
                'kap':{
                    'gkabar':{
                        'origin':'soma'
                    },
                },
                'kdr':{
                    'gkdrbar':{
                        'origin':'soma'
                    },
                },
                'km2':{
                    'gkmbar':{
                        'origin':'soma'
                    },
                },
                'nax_milstein':{
                    'gbar':{
                        'value':0.5*.04
                    },
                    'sh':{
                        'value':1.7
                    },
                },
                'pas':{
                    'e':{
                        'value':-67.
                    },
                    'g':{
                        'origin':'soma'
                    },
                },
            },
            # basal
            #-----------------------------------------
            'basal': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'ions':{
                    'ek':{
                        'origin':'soma'
                    },
                },
                'h_milstein':{
                    'eh':{
                        'origin':'soma'
                    },
                    'ghbar':{
                        'origin':'soma'
                    },
                    'vhalfl':[
                        {'max_loc':75., 'origin':'soma'},
                        {'min_loc':75., 'origin':'soma', 'value':-81}
                    ],
                },
                'kad':{
                    # 'eh':{
                    #     'origin':'soma'
                    # },
                    'gkabar':[
                        {'max_loc':75., 'origin':'soma', 'value':0.},
                        {'min_loc':75., 'max_loc':300., 'origin':'soma', 'slope':.000183, 'value':0.044225},
                        {'min_loc':300., 'origin':'soma', 'value':0.0854},
                    ],
                },
                'kap':{
                    'gkabar':[
                        {'min_loc':75., 'origin':'soma', 'value':0.},
                        {'max_loc':75., 'origin':'soma', 'slope':0.000183},
                    ],
                },
                'kdr':{
                    'gkdrbar':{
                        'origin':'soma'
                    },
                },
                'nas':{
                    'ar':{
                        'value':0.8
                    },
                    'gbar':{
                        'min':0., 'origin':'soma', 'slope':-0.0002,
                    },
                    'sh':{
                        'value':5.
                    },
                },
                'calH':{
                    'gcalbar':{
                            'value':.00125,
                    },
                },
                'pas':{
                    'e':{
                        'value':-67.
                    },
                    'g':{
                        'origin':'soma'
                    },
                },
                'synapse':{
                    'AMPA_KIN':{
                            'gmax':{
                                'value':0.00079
                            }
                    },
                    'NMDA_KIN5':{
                            'gmax':{
                                'value':0.003026
                            }
                    },
                },
            },
            # trunk
            #------------------------------------------------
            'trunk': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'ions':{
                    'ek':{
                        'origin':'soma'
                    },
                },
                'h_milstein':{
                    'eh':{
                        'origin':'soma'
                    },
                    'ghbar':{
                        'origin':'soma',
                        'slope':0.00185,
                        'tau':383.1,
                        'xhalf':260.7,
                    },
                    'vhalfl':[
                        {'max_loc':75., 'origin':'soma'},
                        {'min_loc':75., 'origin':'soma', 'value':-81}
                    ],
                },
                'kad':{
                    # 'eh':{
                    #     'origin':'soma'
                    # },
                    'gkabar':[
                        {'max_loc':75., 'origin':'soma', 'value':0.},
                        {'min_loc':75., 'max_loc':300., 'origin':'soma', 'slope':.000183, 'value':0.044225},
                        {'min_loc':300., 'origin':'soma', 'value':0.0854},
                    ],
                },
                'kap':{
                    'gkabar':[
                        {'min_loc':75., 'origin':'soma', 'value':0.},
                        {'max_loc':75., 'origin':'soma', 'slope':0.000183},
                    ],
                },
                'kdr':{
                    'gkdrbar':{
                        'origin':'soma'
                    },
                },
                'nas':{
                    'ar':{
                        'value':0.8
                    },
                    'gbar':{
                        'origin':'soma',
                    },
                    'sh':{
                        'origin':'soma',
                    },
                    'sha':{
                        'value':5.
                    },
                },
                'pas':{
                    'e':{
                        'origin':'soma',
                        'slope':-1.582714325,
                        'tau':121.9,
                    },
                    'g':{
                        'origin':'soma',
                        'slope':-1.63e-5,
                        'tau':121.9,
                    },
                },
                'synapse':{
                    'AMPA_KIN':{
                            'gmax':{
                                'origin':'soma',
                                'slope':0.000297,
                                'tau':132.8,
                                'value':0.00079,
                            },
                    },
                    'NMDA_KIN5':{
                            'gmax':{
                                'value':0.003026
                            }
                    },
                },
            },
            # apical
            #--------------------------------------------------------
            'apical': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'ions':{
                    'ek':{
                        'origin':'soma'
                    },
                },
                'h_milstein':{
                    'eh':{
                        'origin':'soma'
                    },
                    'ghbar':{
                        'origin':'trunk',
                    },
                    'vhalfl':[
                        {'max_loc':75., 'origin':'soma'},
                        {'min_loc':75., 'origin':'soma', 'value':-81}
                    ],
                },
                'kad':{
                    # 'eh':{
                    #     'origin':'soma'
                    # },
                    'gkabar':[
                        {'max_loc':75., 'origin':'soma', 'value':0.},
                        {'min_loc':75., 'max_loc':300., 'origin':'soma', 'slope':.000183, 'value':0.044225},
                        {'min_loc':300., 'origin':'soma', 'value':0.0854},
                    ],
                },
                'kap':{
                    'gkabar':[
                        {'min_loc':75., 'origin':'soma', 'value':0.},
                        {'max_loc':75., 'origin':'soma', 'slope':0.000183},
                    ],
                },
                'kdr':{
                    'gkdrbar':{
                        'origin':'soma'
                    },
                },
                'nas':{
                    'ar':{
                        'value':0.8
                    },
                    'gbar':{
                        'origin':'trunk',
                        'min':0.0,
                        'slope':-0.0002,
                    },
                    'sh':{
                        'origin':'soma',
                    },
                    'sha':{
                        'value':5.
                    },
                },
                'pas':{
                    'e':{
                        'origin':'trunk',
                    },
                    'g':{
                        'origin':'trunk',
                    },
                },
                'synapse':{
                    'AMPA_KIN':{
                            'gmax':{
                                'origin':'trunk',
                            },
                    },
                    'NMDA_KIN5':{
                            'gmax':{
                                'origin':'trunk',
                            }
                    },
                },
            },
            # tuft
            #---------------------------------------------
            'tuft': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'ions':{
                    'ek':{
                        'origin':'soma'
                    },
                },
                'h_milstein':{
                    'eh':{
                        'origin':'soma'
                    },
                    'ghbar':{
                        'origin':'trunk',
                    },
                    'vhalfl':[
                        {'max_loc':75., 'origin':'soma'},
                        {'min_loc':75., 'origin':'soma', 'value':-81}
                    ],
                },
                'kad':{
                    # 'eh':{
                    #     'origin':'soma'
                    # },
                    'gkabar':[
                        {'max_loc':75., 'origin':'soma', 'value':0.},
                        {'min_loc':75., 'max_loc':300., 'origin':'soma', 'slope':.000183, 'value':0.044225},
                        {'min_loc':300., 'origin':'soma', 'value':0.0854},
                    ],
                },
                'kap':{
                    'gkabar':[
                        {'min_loc':75., 'origin':'soma', 'value':0.},
                        {'max_loc':75., 'origin':'soma', 'slope':0.000183},
                    ],
                },
                'kdr':{
                    'gkdrbar':{
                        'origin':'soma'
                    },
                },
                'nas':{
                    'ar':{
                        'value':0.8
                    },
                    'gbar':{
                        'origin':'trunk',
                    },
                    'sh':{
                        'origin':'soma',
                    },
                    'sha':{
                        'value':5.
                    },
                },
                'pas':{
                    'e':{
                        'origin':'trunk',
                    },
                    'g':{
                        'origin':'trunk',
                    },
                },
                'synapse':{
                    'AMPA_KIN':{
                            'gmax':{
                                'origin':'trunk',
                            },
                    },
                    'NMDA_KIN5':{
                            'gmax':{
                                'origin':'trunk',
                            }
                    },
                },
            },
            # spine_head
            #-------------------------------------------
            'spine_head': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'pas':{
                    'e':{
                        'value':-67,
                    },
                    'g':{
                        'origin':'parent',
                    },
                },
            },
            # spine_neck
            #----------------------------------------------
            'spine_neck': {
                'cable':{
                    'Ra':{
                        'origin':'soma'
                    },
                    'cm':{
                        'origin':'soma'
                    },
                },
                'pas':{
                    'e':{
                        'value':-67,
                    },
                    'g':{
                        'origin':'parent',
                    },
                },
            },
}