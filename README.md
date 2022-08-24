# Realistic-Image-Synthesis
Here, in this course we learned the advanced concepts of computer graphics 1 at Saarland University. The concepts we learned are:<br/>
1. Monte Carlo Sampling <br/>
2. Vanilla Path Tracing <br/>
| Glossy Surface (100 spp)     | Diffuse Surface (100spp) |
| ----------- | ----------- |
| <img src="imgs/PT/cornell_box_glossy_100.png" width="300"/>      | <img src="imgs/PT/cornell_box_100.png" width="300"/>       |
3. Next Event Estimation <br/>
Combination of direect lighting and path tracing.
Direct lighting:Connect to a random light sources at each hit point and collect the contributions as the ray bounces.<br/> 
Path Tracing : Collect the throughput as ray bounces off. 
MIS(Multiple Importance sampling): Combine the NEE with path tracing more effectively as NEE converge faster for some lighting conditions.
|Surfaces |NEE     | MIS |
|----------- |----------- | ----------- |
| Glossy(10 spp)|<img src="imgs/NEE/NEE/NEE_cornel_box_glossy_10spp.png" width="300"/>      | <img src="imgs/NEE/MIS/MIS_cornel_box_glossy_10spp.png" width="300"/>       |
| Diffuse(10 spp)|<img src="imgs/NEE/NEE/NEE_cornel_box_10spp.png" width="300"/>      | <img src="imgs/NEE/MIS/MIS_cornel_box_10spp.png" width="300"/>       |
| Water(200 spp)|<img src="imgs/NEE/NEE/NEE_cornel_box_water.png" width="300"/>      | <img src="imgs/NEE/MIS/MIS_cornel_box_water_200spp.png" width="300"/>       |

6. Photon Mapping and Density Estimation <br/>
