# Realistic-Image-Synthesis
Here, in this course we learned the advanced concepts of computer graphics 1 at Saarland University. The concepts we learned are:<br/>
## 1. Monte Carlo Sampling <br/>
## 2. Vanilla Path Tracing <br/>

| Glossy Surface (100 spp) | Diffuse Surface (100spp) |
| :-----------: | :-----------: |
 <img src="imgs/PT/cornell_box_glossy_100.png" width="300"/>      | <img src="imgs/PT/cornell_box_100.png" width="300"/>  

## 3. Next Event Estimation <br/>
Combination of **direect lighting** and **path tracing**.<br/>
**Direct lighting**:Connect to a random light sources at each hit point and collect the contributions as the ray bounces.<br/> 
**Path Tracing** : Collect the throughput as ray bounces off. <br/>
MIS(Multiple Importance sampling): Combine the NEE with path tracing more effectively as NEE converge faster for some lighting conditions.<br/>

Surfaces | NEE | MIS
:-----: | :----: | :-----:
Glossy(10 spp)   | <img src="imgs/NEE/NEE/NEE_cornel_box_glossy_10spp.png" width="300"/> | <img src="imgs/NEE/MIS/MIS_cornel_box_glossy_10spp.png" width="300"/>
Diffuse(10 spp)  | <img src="imgs/NEE/NEE/NEE_cornel_box_10spp.png" width="300"/>  | <img src="imgs/NEE/MIS/MIS_cornel_box_10spp.png" width="300"/>
Water(200 spp)   | <img src="imgs/NEE/NEE/NEE_cornel_box_water.png" width="300"/>  | <img src="imgs/NEE/MIS/MIS_cornel_box_water_200spp.png" width="300"/>


## 4. Photon Mapping and Density Estimation <br/>
Combination of **photon tracing** and **Density estimation**.<br/>
**Photon Tracing** : Shoot photons from light sources and store them at the hit pts using some acceleration structure for computing the density during density estimation. Remember to not store photons on specular surfaces. <br/>
**Density Estimation** : Estimate the density around a point taking the contribution of all the photons in the neighbourhood of this point using some kernel (we have used epanechnikov kernel) with certain radius. We didn't perform density estimation on specular surfaces as the probablility of sampling a direction in the exact direction (in the specular lobe) is very low.<br/> 
**Direct Lighting with PM** : Glossy materials don't work well. Hence, to improve efficiency of progressive photon mapping in the presence of glossy materials is
to not do density estimation on glossy materials either. Instead, compute the direct illumination at a glossy hit point, and bounce, as you would do in path tracing. Keep in mind that you will have to either ignore hitting the light source after a glossy bounce,
or use MIS to combine direct illumination with randomly hitting the light. <br/>

Surfaces | PM | MIS
:-----: | :----: | :-----:
Water(100 spp)   | <img src="imgs/PM/Density_estimation_7.3(fig1)/cornel_box_water_100spp.png" width="300"/>  | <img src="imgs/PM/with_glossy_improved_7.5(fig3)/cornel_box_water_improved_100spp.png" width="300"/>
Specular(100 spp)  | <img src="imgs/PM/Density_estimation_7.3(fig1)/cornel_box_specular_100spp.png" width="300"/>  | same as no glossy surfaces

## 5. HDR and Tone Mapping <br/>
**HDR(High Dynamic Range)**: It captures all the dynamic range of an image.<br/>
**TM(Tone Mapping)**: Convert the high dynamic range to displayable range without losing any details as many viewing devices are not able to display all the ranges.<br/>
Here, we have computed HDR of two images using multi image exposure technique and performed tone mapping using 4 tone mapping techniques.<br/>
Operators | Day Image | Night Image
:-----: | :----: | :-----:
Drago   | <img src="imgs/HDR_and_TM/daytime/drago/TMO_DragoTMO_dLdMax120_db0.65.png" width="300"/>  | <img src="imgs/HDR_and_TM/nighttime/drago/TMO_DragoTMO_dLdMax80_db0.65.png" width="300"/>
Durand  | <img src="imgs/HDR_and_TM/nighttime/durand/TMO_Durand_contrast7.png" width="300"/>  | <img src="imgs/HDR_and_TM/nighttime/durand/TMO_Durand_contrast7.png" width="300"/>
Logarithmic   | <img src="imgs/HDR_and_TM/daytime/logarithmic/TMO_Logarithmic_logQ0.5_logK0.5.png" width="300"/>  | <img src="imgs/HDR_and_TM/nighttime/logarithmic/TMO_Logarithmic_logQ0.5_logK0.5.png" width="300"/>
Reinhard  | <img src="imgs/HDR_and_TM/daytime/reinhard/Bottles_Small_TMO_Reinhard_palpha0.18_pwhite2.5_local.png" width="300"/>  | <img src="imgs/HDR_and_TM/nighttime/reinhard/Bottles_Small_TMO_Reinhard_palpha0.18_pwhite2.5_local.png" width="300"/>
