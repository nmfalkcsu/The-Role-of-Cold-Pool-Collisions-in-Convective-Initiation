module sounding_SF12
use constants
    use constants_b

    implicit none

    

    real, dimension(NZ) :: &
        THETAB_Z = (/301.1853,&
301.1853,&
301.1894,&
301.195,&
301.225,&
301.2128,&
301.2196,&
301.2856,&
301.3362,&
301.4428,&
301.5674,&
301.6487,&
301.7577,&
301.9127,&
302.0226,&
302.0911,&
302.2322,&
302.428,&
302.7463,&
302.9108,&
303.1788,&
303.404,&
303.8825,&
304.1006,&
304.3102,&
304.3672,&
304.4206,&
304.7224,&
305.0233,&
305.1503,&
305.269,&
305.5364,&
305.7685,&
305.9409,&
306.1588,&
306.4125,&
306.5974,&
306.8407,&
307.1887,&
307.4681,&
307.6466,&
307.9165,&
308.1319,&
308.5438,&
309.1626,&
309.5025,&
309.9938,&
311.2311,&
311.5846,&
311.8589,&
312.395,&
312.564,&
312.8203,&
313.1745,&
313.4961,&
313.8107,&
314.113,&
314.2235,&
314.2601,&
314.8564,&
315.1808,&
315.43,&
315.5689,&
315.8422,&
316.1181,&
316.5173,&
316.7131,&
316.8874,&
317.1982,&
317.7656,&
318.2406,&
318.4812,&
318.7886,&
319.2188,&
319.3408,&
319.4971,&
319.817,&
320.0953,&
320.2407,&
320.4636,&
320.5431,&
320.5717,&
320.5675,&
320.732,&
320.9449,&
321.1008,&
321.2966,&
321.5022,&
321.6204,&
321.6469,&
321.7955,&
321.9068,&
322.0641,&
322.1526,&
322.2548,&
322.3363,&
322.3999,&
322.6875,&
322.8883,&
323.2643,&
323.5302,&
323.6847,&
323.6937,&
323.7688,&
323.9603,&
324.1159,&
324.3338,&
324.5136,&
324.7035,&
325.0412,&
325.2402,&
325.3795,&
325.488,&
325.6698,&
325.8092,&
325.9074,&
325.9975,&
326.1045,&
326.2102,&
326.4161,&
326.6984,&
327.0252,&
327.3882,&
327.7512,&
328.1142,&
328.4772,&
328.8402,&
329.2032,&
329.5662,&
329.9292,&
330.2922,&
330.6552,&
331.0183,&
331.3813,&
331.7443,&
332.0866,&
332.3208,&
332.5549,&
332.7891,&
333.0232,&
333.2574,&
333.4915,&
333.7257,&
333.9598,&
334.194,&
334.4281,&
334.6623,&
334.8964,&
335.1306,&
335.3647,&
335.5989,&
335.833,&
336.0672,&
336.2953,&
336.459,&
336.6228,&
336.7865,&
336.9503,&
337.1141,&
337.2778/),&

rvb_z = (/0.019523,&
0.019523,&
0.019465,&
0.019487,&
0.019251,&
0.019325,&
0.01932,&
0.019111,&
0.018848,&
0.018553,&
0.018109,&
0.017806,&
0.017634,&
0.017153,&
0.017117,&
0.017099,&
0.016888,&
0.016578,&
0.015852,&
0.015579,&
0.01502,&
0.01515,&
0.014507,&
0.01431,&
0.013722,&
0.013639,&
0.01413,&
0.013394,&
0.013126,&
0.013489,&
0.013337,&
0.012948,&
0.012656,&
0.012674,&
0.012198,&
0.012618,&
0.012296,&
0.012235,&
0.011809,&
0.010983,&
0.011115,&
0.011009,&
0.010334,&
0.011085,&
0.0098071,&
0.0089784,&
0.0090022,&
0.0044701,&
0.0055321,&
0.0078202,&
0.00636,&
0.0064638,&
0.0061335,&
0.0056101,&
0.0056781,&
0.0056077,&
0.0057107,&
0.0057913,&
0.0076351,&
0.0077336,&
0.0074274,&
0.0074444,&
0.007536,&
0.0066616,&
0.0064143,&
0.0057791,&
0.0066208,&
0.0070249,&
0.0070008,&
0.0068286,&
0.0064562,&
0.0065491,&
0.0064856,&
0.0055554,&
0.0051994,&
0.0050598,&
0.0047037,&
0.0045857,&
0.0045622,&
0.004313,&
0.0041681,&
0.0041234,&
0.0040917,&
0.0037769,&
0.0035788,&
0.0037549,&
0.0032165,&
0.0028402,&
0.002551,&
0.0026277,&
0.0025978,&
0.0026662,&
0.0030188,&
0.00316,&
0.0030592,&
0.002691,&
0.0030779,&
0.0027929,&
0.0029337,&
0.0025851,&
0.0019565,&
0.0020866,&
0.0023256,&
0.0023241,&
0.0025128,&
0.0026031,&
0.0022576,&
0.0018872,&
0.0018799,&
0.0017907,&
0.0018294,&
0.0020024,&
0.0021999,&
0.0022282,&
0.002375,&
0.0025112,&
0.0025561,&
0.0027498,&
0.0028443,&
0.002857,&
0.0028493,&
0.0028356,&
0.0028385,&
0.0028414,&
0.0028442,&
0.0028471,&
0.00285,&
0.0028529,&
0.0028558,&
0.0028587,&
0.0028616,&
0.0028645,&
0.0028674,&
0.0028703,&
0.0028732,&
0.0028695,&
0.0028318,&
0.002794,&
0.0027563,&
0.0027185,&
0.0026808,&
0.002643,&
0.0026053,&
0.0025675,&
0.0025297,&
0.002492,&
0.0024542,&
0.0024165,&
0.0023787,&
0.002341,&
0.0023032,&
0.0022655,&
0.0022277,&
0.0021901,&
0.0021536,&
0.0021171,&
0.0020806,&
0.0020441,&
0.0020077,&
0.0019712/),&

thetab_w_z = (/301.1158,&
301.1158,&
301.1909,&
301.1973,&
301.221,&
301.2127,&
301.2134,&
301.2269,&
301.3253,&
301.3777,&
301.5078,&
301.597,&
301.6735,&
301.8535,&
301.9802,&
302.0534,&
302.1571,&
302.3192,&
302.5657,&
302.8398,&
303.0199,&
303.2631,&
303.6439,&
303.9944,&
304.2348,&
304.3286,&
304.3161,&
304.6233,&
304.8295,&
305.0784,&
305.2237,&
305.4021,&
305.6401,&
305.8933,&
305.9902,&
306.2743,&
306.5183,&
306.6821,&
307.0671,&
307.2983,&
307.5681,&
307.8148,&
307.9833,&
308.2447,&
308.9095,&
309.3596,&
309.6835,&
310.7082,&
311.4525,&
311.623,&
312.2412,&
312.4631,&
312.7077,&
312.9979,&
313.3291,&
313.6531,&
313.9912,&
314.1867,&
314.2118,&
314.5237,&
315.0475,&
315.2827,&
315.4997,&
315.7145,&
315.9284,&
316.3485,&
316.6195,&
316.7963,&
317.0382,&
317.4584,&
318.0734,&
318.3392,&
318.6173,&
319.0379,&
319.3061,&
319.4074,&
319.6401,&
319.9911,&
320.1512,&
320.3528,&
320.5265,&
320.556,&
320.5814,&
320.5962,&
320.8633,&
321.0127,&
321.2023,&
321.3684,&
321.6051,&
321.6236,&
321.7056,&
321.8462,&
321.9775,&
322.1092,&
322.2106,&
322.3021,&
322.3518,&
322.5255,&
322.8005,&
323.0581,&
323.4128,&
323.6495,&
323.6827,&
323.7122,&
323.8355,&
324.0549,&
324.2152,&
324.4412,&
324.5992,&
324.8662,&
325.1462,&
325.3127,&
325.4327,&
325.5733,&
325.7376,&
325.8635,&
325.9527,&
326.0505,&
326.1632,&
326.2703,&
326.5573,&
326.8437,&
327.2067,&
327.5697,&
327.9327,&
328.2957,&
328.6587,&
329.0217,&
329.3847,&
329.7477,&
330.1107,&
330.4737,&
330.8368,&
331.1998,&
331.5628,&
331.9258,&
332.2037,&
332.4378,&
332.672,&
332.9061,&
333.1403,&
333.3744,&
333.6086,&
333.8427,&
334.0769,&
334.311,&
334.5452,&
334.7794,&
335.0135,&
335.2477,&
335.4818,&
335.716,&
335.9501,&
336.1843,&
336.3772,&
336.5409,&
336.7047,&
336.8684,&
337.0322,&
337.1959/),&

rvb_w_z = (/0.019573,&
0.019573,&
0.01948,&
0.019387,&
0.019277,&
0.019278,&
0.019321,&
0.019322,&
0.018876,&
0.018657,&
0.018344,&
0.018035,&
0.017812,&
0.017292,&
0.017056,&
0.017156,&
0.017033,&
0.016709,&
0.016345,&
0.015645,&
0.0155,&
0.01484,&
0.014957,&
0.014361,&
0.01381,&
0.013841,&
0.014348,&
0.013456,&
0.013665,&
0.013337,&
0.013391,&
0.013355,&
0.012629,&
0.012796,&
0.012589,&
0.012606,&
0.012311,&
0.012241,&
0.011962,&
0.011621,&
0.010629,&
0.011078,&
0.010943,&
0.010796,&
0.010368,&
0.0092827,&
0.0090435,&
0.00671,&
0.0047817,&
0.0075006,&
0.0069176,&
0.0062787,&
0.0064012,&
0.0058141,&
0.0056696,&
0.005731,&
0.005398,&
0.0059824,&
0.0063306,&
0.0076132,&
0.0075535,&
0.0074983,&
0.0073986,&
0.0072497,&
0.0066394,&
0.0057908,&
0.0060541,&
0.0070164,&
0.0066973,&
0.0070059,&
0.0064414,&
0.0065048,&
0.0066163,&
0.0060552,&
0.0051706,&
0.0050907,&
0.0049352,&
0.0045855,&
0.0045286,&
0.0044269,&
0.0041852,&
0.0041441,&
0.0041057,&
0.004096,&
0.0033944,&
0.0036182,&
0.003587,&
0.0030509,&
0.0025223,&
0.0025741,&
0.0025862,&
0.0026668,&
0.0029556,&
0.0031029,&
0.0031743,&
0.0028108,&
0.0027605,&
0.0032091,&
0.0026022,&
0.0029255,&
0.0021828,&
0.001864,&
0.0023086,&
0.0023558,&
0.0023483,&
0.0026547,&
0.0024656,&
0.0019777,&
0.0018648,&
0.0019009,&
0.0017822,&
0.0019025,&
0.0020972,&
0.0022607,&
0.0022702,&
0.002483,&
0.0025388,&
0.0026044,&
0.002782,&
0.0028627,&
0.0028598,&
0.0028341,&
0.002837,&
0.0028399,&
0.0028428,&
0.0028457,&
0.0028486,&
0.0028515,&
0.0028544,&
0.0028573,&
0.0028602,&
0.0028631,&
0.0028659,&
0.0028688,&
0.0028717,&
0.0028746,&
0.0028507,&
0.0028129,&
0.0027752,&
0.0027374,&
0.0026996,&
0.0026619,&
0.0026241,&
0.0025864,&
0.0025486,&
0.0025109,&
0.0024731,&
0.0024354,&
0.0023976,&
0.0023599,&
0.0023221,&
0.0022843,&
0.0022466,&
0.0022088,&
0.0021718,&
0.0021353,&
0.0020989,&
0.0020624,&
0.0020259,&
0.0019894/)

end module sounding_SF12