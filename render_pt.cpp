#include "../scene.h"
#include "../color.h"
#include "../samplers.h"
#include "../cameras.h"
#include "../hash.h"
#include "../debug.h"
#include "../renderer.h"

/// Path Tracing with MIS and Russian Roulette.
class PathTracingRenderer : public Renderer {
public:
    PathTracingRenderer(const Scene& scene, size_t max_path_len)
        : Renderer(scene), max_path_len(max_path_len)
    {}

    std::string name() const override { return "pt"; }

    void reset() override { iter = 1; }

    void render(Image& img) override {
        auto kx = 2.0f / (img.width - 1);
        auto ky = 2.0f / (img.height - 1);

        process_tiles(0, 0, img.width, img.height,
            default_tile_width, default_tile_height,
            [&](size_t xmin, size_t ymin, size_t xmax, size_t ymax) {
                UniformSampler sampler(sampler_seed(xmin ^ ymin, iter));
                for (size_t y = ymin; y < ymax; y++) {
                    for (size_t x = xmin; x < xmax; x++) {
                        auto ray = scene.camera->gen_ray(
                            (x + sampler()) * kx - 1.0f,
                            1.0f - (y + sampler()) * ky);

                        debug_raster(x, y);
                        img(x, y) += rgba(path_trace(ray, sampler), 1.0f);
                    }
                }
            });
        iter++;
    }

    inline rgb path_trace(Ray ray, Sampler& sampler);

private:
    size_t max_path_len;
    size_t iter;
};

rgb PathTracingRenderer::path_trace(Ray ray, Sampler& sampler) {
    rgb color(0.0f);
    float surv_prob = 0.0f;
    rgb beta(1.0f);

    auto dir_lit = false; //specular
    float weight_brdf = 0.0f;
    float3 beta_brdf = float3(1.0f);
    float pdf_nee = 0.0f;
    BsdfSample bsdf_new;
    //float pdf_brdf = 1.0f;
    for (size_t path_len = 0; path_len < max_path_len; path_len++) {
        Hit hit = scene.intersect(ray);
        if (hit.tri < 0) break;

        auto surf = scene.surface_params(ray, hit);
        auto mat = scene.material(hit);
        auto out = -ray.dir;
        
        if (auto light = mat.emitter) {
            // Direct hits on a light source
            if (surf.entering) {
                // TODO: compute the incoming radiance from the emitter.           
                auto light_head = light->emission(out, sampler(), sampler());


                if (dir_lit ||path_len == 0) {
                    color = color + (beta)*light_head.intensity;
                }
    
                // MIS brdf
                if (path_len != 0) {
                    auto nee_brdf_length = length(surf.point - ray.org);
                    auto nee_brdf_cos = dot((-surf.point + ray.org) / nee_brdf_length, surf.face_normal);
                    auto pdf_nee_brdf = light_head.pdf_area * (nee_brdf_length * nee_brdf_length) / nee_brdf_cos;

                    auto pdf_brdf = bsdf_new.pdf;
                    auto weight_brdf = pdf_brdf / (pdf_brdf + pdf_nee_brdf);
                    auto beta_brdf = bsdf_new.color / pdf_brdf;
                    color = color + (beta) * weight_brdf * beta_brdf * light_head.intensity;
                }
               
                return color;
            }               
        }

        // Materials without BSDFs act like black bodies
        if (!mat.bsdf) break;

        bool specular = mat.bsdf->type() == Bsdf::Type::Specular;

        // TODO: Evaluate direct lighting
        // When using Next Event Estimation, you should select a light in the scene
        // (uniformly, for now) and compute direct illumination from it.
        // Important notes:
        //   - Weight the contribution of one light with the probability of choosing that light,
        //   - Be careful when combining Next Event Estimation with BRDF sampling (see assignment),
        //   - Do not take lights that face away the surface at the hit point into account.

        //Sample new direction
        bsdf_new = mat.bsdf->sample(sampler, surf, out,specular);
        
        //direct illumination

        // select a random light index
        auto random_lt_no = int(sampler() * scene.lights.size());
        // make sure it lies with in the array lenth of lights
        random_lt_no = std::min(random_lt_no, int(scene.lights.size() - 1));
        // Apply sample direct on the light whose index we calculated
        auto sample_dir_lt = scene.lights[random_lt_no]->sample_direct(surf.point, sampler);
        auto light_hs_area = scene.lights[random_lt_no]->has_area();
        // Compute shadow ray from aurface hitpoint to the selected light pos
        auto dist = length(sample_dir_lt.pos - surf.point);
        auto shadow_ray = Ray(surf.point, (sample_dir_lt.pos - surf.point), offset, 1.0f - offset);
        

        //cosine of the surface normal and the light direction or visibility
        auto cos_x = dot(shadow_ray.dir / dist, surf.face_normal);
        

        if (!scene.occluded(shadow_ray)) {
            if (!specular) {
                // comput pdf for nee for different lights (area,point)
                if (light_hs_area) {
                    pdf_nee = (sample_dir_lt.pdf_area * (dist * dist)) / sample_dir_lt.cos;

                }
                else {
                    pdf_nee = (dist * dist) / sample_dir_lt.cos;

                }

                // for MIS of the NEE
                auto weight_nee = pdf_nee / (pdf_nee + mat.bsdf->pdf(shadow_ray.dir / dist, surf, ray.dir));

                auto beta_nee = (mat.bsdf->eval(ray.dir, surf, shadow_ray.dir / dist) * scene.lights.size() * cos_x) / (pdf_nee);

                if (cos_x > 0.0f) {
                    color = color + beta * beta_nee * weight_nee * sample_dir_lt.intensity;
                }
            }

        }

        dir_lit = specular;
        //break; // uncomment this debugging direct lighting

        // TODO: Add Russian Roulette to prevent infinite recursion        
        if (path_len > 3) {
            auto surv_prob = russian_roulette(beta, 0.75f);
            if (sampler() <= surv_prob) {
                beta = beta / surv_prob;
            }
            else {
                break;
            }
        }


        //path contribution
        beta *= (bsdf_new.color) / (bsdf_new.pdf);

        //Trace new ray
        ray = Ray(surf.point, normalize(bsdf_new.in), offset);
        

         //if ray don't hit light and reach max path length then send black color.
        if (path_len == max_path_len-1) {
            return rgb(0.0f);
        }

        // TODO:
        // If the path is not terminated, sample a direction to continue the path with from the material.
        // You can do this using the sample() function of the Bsdf class. Update the path weight, and do not
        // forget to take the Russian Roulette probability into account!


        // General remarks:
        //   - The algorithm is currently expressed in an iterative form.
        //     If you do not feel confident with it, feel free to use a recursive form instead.
        //   - Because this function is called in parallel, make sure you avoid data races.
        //     Disable OpenMP if you are encountering any issue.
        //   - The emission() function in the Light class is required to evaluate the radiance when hitting a light source.
        //     Of all the fields in the EmissionValue structure, you will only need two: pdf_area (only for MIS) and intensity.
        //   - The sample_direct() function in the Light class is required to evaluate direct illumination at a given point.
        //     You will only need the pdf_area (only for MIS), cos, pos, and intensity fields.
    }
    return color;
}

std::unique_ptr<Renderer> create_pt_renderer(const Scene& scene, size_t max_path_len) {
    return std::unique_ptr<Renderer>(new PathTracingRenderer(scene, max_path_len));
}
