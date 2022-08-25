#include "../scene.h"
#include "../color.h"
#include "../samplers.h"
#include "../cameras.h"
#include "../hash.h"
#include "../hash_grid.h"
#include "../debug.h"
#include "../renderer.h"

#ifdef USE_STD_THREAD
#include <mutex>
#endif

#ifdef _OPENMP
#include <omp.h>
#define OMP_THREAD_NUM() omp_get_thread_num()
#else
#define OMP_THREAD_NUM() 0
#endif

struct Photon {
    rgb contrib;         ///< Path contribution
    SurfaceParams surf;  ///< Surface parameters at the vertex
    float3 in_dir;       ///< Incoming direction

    Photon() {}
    Photon(const rgb& c, const SurfaceParams& s, const float3& i)
        : contrib(c), surf(s), in_dir(i)
    {}
};

class PhotonMappingRenderer : public Renderer {
public:
    PhotonMappingRenderer(const Scene& scene, size_t max_path_len)
        : Renderer(scene), max_path_len(max_path_len)
    {}

    std::string name() const override { return "ppm"; }

    void reset() override { iter = 1; }

    void render(Image& img) override {
        constexpr float alpha = 0.75f;

        static float base_radius = 1.0f;
        if (iter == 1)
            base_radius = 2.0f * estimate_pixel_size(img.width, img.height);

        auto kx = 2.0f / (img.width - 1);
        auto ky = 2.0f / (img.height - 1);

        auto light_path_count = img.width * img.height;
        photons.clear();

        // Trace a light path for every pixel
#ifdef USE_STD_THREAD
        std::mutex photon_mutex;
        size_t batch_size = 32;
        parallel_for(0, light_path_count / batch_size + (light_path_count % batch_size ? 1 : 0),
        [&](size_t batch_id) {
            std::vector<Photon> photon_buffer;
            UniformSampler sampler(sampler_seed(batch_id, iter));

            int num = std::min((batch_id + 1) * batch_size, light_path_count) - batch_id * batch_size;
            for (int i = 0; i < num; ++i)
                trace_photons(photon_buffer, sampler);

            { std::lock_guard guard(photon_mutex);
                photons.insert(photons.end(), photon_buffer.begin(), photon_buffer.end());
            }
        });
#else
        #pragma omp parallel
        {
            std::vector<Photon> photon_buffer;
            UniformSampler sampler(sampler_seed(OMP_THREAD_NUM(), iter));

            #pragma omp for schedule(dynamic) nowait
            for (size_t i = 0; i < light_path_count; i++)
                trace_photons(photon_buffer, sampler);

            #pragma omp critical
            { photons.insert(photons.end(), photon_buffer.begin(), photon_buffer.end()); }
        }
#endif

        // Build the photon map
        radius = base_radius / std::pow(float(iter), 0.5f * (1.0f - alpha));
        photon_map.build(
            [&] (size_t i) { return photons[i].surf.point; },
            photons.size(),
            radius);

        // Trace the eye paths
        process_tiles(0, 0, img.width, img.height,
            default_tile_width, default_tile_height,
            [&] (size_t xmin, size_t ymin, size_t xmax, size_t ymax) {
            UniformSampler sampler(sampler_seed(xmin ^ ymin, iter));
            for (size_t y = ymin; y < ymax; y++) {
                for (size_t x = xmin; x < xmax; x++) {
                    auto ray = scene.camera->gen_ray((x + sampler()) * kx - 1.0f, 1.0f - (y + sampler()) * ky);
                    debug_raster(x, y);
                    img(x, y) += atomically(rgba(trace_eye_path(ray, sampler, light_path_count), 1.0f));
                }
            }
        });
        iter++;
    }

    void trace_photons(std::vector<Photon>& photons, Sampler& sampler);
    rgb trace_eye_path(Ray ray, Sampler& sampler, size_t light_path_count);
    float estimate_pixel_size(size_t w, size_t h);

private:
    std::vector<Photon> photons;
    HashGrid photon_map;
    size_t max_path_len;
    size_t iter;
    float radius;
};

void PhotonMappingRenderer::trace_photons(std::vector<Photon>& photons, Sampler& sampler) {
    // TODO: Choose a light to sample from (uniformly) and get an emission sample for it
    // Hints:
    // _ Lights are in scene.lights
    // _ Emission samples can be obtained by calling sample_emission() on a given light

    auto random_lt_no = int(sampler() * scene.lights.size());
    random_lt_no = std::min(random_lt_no, int(scene.lights.size() - 1));
    auto sample_dir_lt = scene.lights[random_lt_no]->sample_emission(sampler);

    
    // TODO: Create the starting ray from the light sample
    
    // Hints:
    // _ Add an offset to avoid artifacts. The constructor for Ray is: Ray(origin, direction, offset)
    // _ Initialize the pdf of the path according to the emission sample
    // _ Initialize the contribution of the path according to the emission sample
    auto ray = Ray(sample_dir_lt.pos, sample_dir_lt.dir, offset); // TODO: initialize this correctly!
    //pdf of path and pdf of light as we are sampling light as well as path.
    float pdf_path = sample_dir_lt.pdf_area * sample_dir_lt.pdf_dir ; 
    auto  beta = sample_dir_lt.intensity * scene.lights.size() * sample_dir_lt.cos / pdf_path;  //path contribution

    for (size_t path_len = 0; path_len < max_path_len; path_len++) {
        Hit hit = scene.intersect(ray);
        if (hit.tri < 0) break;

        auto mat = scene.material(hit);
        auto surf = scene.surface_params(ray, hit);
        auto out = -ray.dir;
        if (!mat.bsdf) break;

        // TODO: Implement photon shooting here
        // Hints:
        // _ Do not store photons on specular surfaces (bsdf->type() == Bsdf::Type::Specular)
        // _ Use Russian Roulette to terminate photon paths (take it into account in the pdf)
        // _ Add photons to the "photons" vector
        /*break;*/ // TODO: Remove 
        bool specular = mat.bsdf->type() == Bsdf::Type::Specular;

        if (!specular) {
            photons.push_back(Photon(beta, surf, out));
        }
        // For the new direction of photon
        auto bsdf_new = mat.bsdf->sample(sampler, surf, out, specular);


        //Russian roullette

        auto surv_prob = russian_roulette(beta, 0.75f);
        if (sampler() <= surv_prob) {
            beta = beta / surv_prob;
        }
        else {
            break;
        }

        //path contribution
        beta *= bsdf_new.color / (bsdf_new.pdf);

        //Trace photon in new path               
        ray = Ray(surf.point, bsdf_new.in, offset);
    }
}

rgb PhotonMappingRenderer::trace_eye_path(Ray ray, Sampler& sampler, size_t light_path_count) {
    static size_t max_path_len = 10;

    // TODO: Initialize path variables (see Path Tracing assignment)
    rgb color(0.0f); // TODO: initialize correctly
    rgb beta_eye(1.0f);
    float pdf_nee = 1.0f;
    auto dir_lit = Bsdf::Type::Diffuse; //specular

    ray.tmin = offset;
    for (size_t path_len = 0; path_len < max_path_len; path_len++) {
        Hit hit = scene.intersect(ray);
        if (hit.tri < 0) break;

        auto surf = scene.surface_params(ray, hit);
        auto mat = scene.material(hit);
        auto out = -ray.dir;

        // TODO: Handle direct light hits (see Path Tracing assignment)
        if (auto light = mat.emitter) {
            // Direct hits on a light source
            if (surf.entering) {
                // TODO: compute the incoming radiance from the emitter.           
                auto light_head = light->emission(out, surf.uv.x, surf.uv.y);
                //Prevent light hit in case of GLossy
                if (dir_lit != Bsdf::Type::Glossy) {
                    color = color + (beta_eye)*light_head.intensity;
                }

                return color;
            }
        }
        if (!mat.bsdf) break;

        // TODO: Do a photon query if the material is not specular, otherwise bounce (as in Path Tracing)
        // Hints:
        // _ Photons can be queried from the "photon_map" object. Use the following function call:
        bool specular = mat.bsdf->type() == Bsdf::Type::Specular;
        //Sample new direction
        auto bsdf_new = mat.bsdf->sample(sampler, surf, out, specular);
        float3 point = surf.point; // TODO: initialize correctly
        if (mat.bsdf->type() == Bsdf::Type::Diffuse) {
            dir_lit = Bsdf::Type::Diffuse;
            photon_map.query(point,
                [&](size_t i) { return photons[i].surf.point; },
                [&](size_t i, float d2) {
                    const Photon& p = photons[i];
                    
                        // p is the photon, d2 is the squared distance from the point to the photon
                        // TODO: Do something
                    
                    auto brdf_term = mat.bsdf->eval(out, surf, p.in_dir);
                    auto k_inner = (1.0f - (d2 / (radius * radius))) * p.contrib;
                    auto kernel_term = 2.0f * k_inner / (pi * radius * radius);
                    color += kernel_term * brdf_term * beta_eye / photons.size();
                });
            return color;
            // _ Use the Epanechnikov filter as the density estimation kernel  
        }

        else if (mat.bsdf->type() == Bsdf::Type::Glossy) {
            dir_lit = Bsdf::Type::Glossy;
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

            pdf_nee = (sample_dir_lt.pdf_area * (dist * dist)) / sample_dir_lt.cos;

            if (!scene.occluded(shadow_ray)) {

                auto beta_nee = (mat.bsdf->eval(ray.dir, surf, shadow_ray.dir / dist) * scene.lights.size() * cos_x) / (pdf_nee);

                if (cos_x > 0.0f) {
                    color = color + beta_eye * beta_nee * sample_dir_lt.intensity;
                }

            }
            // Trace new path and update path weigth
            auto bsdf_sample = mat.bsdf->sample(sampler, surf, out, specular);
            beta_eye *= (bsdf_sample.color) / (bsdf_sample.pdf);
            ray = Ray(surf.point, bsdf_sample.in, offset);

        }
 
        else {
            // If specular
            dir_lit = Bsdf::Type::Specular;
            //path contribution
            beta_eye *= (bsdf_new.color) / (bsdf_new.pdf);
            //Trace new ray
            ray = Ray(surf.point, bsdf_new.in, offset);
        }


    }

    return color;
}

float PhotonMappingRenderer::estimate_pixel_size(size_t w, size_t h) {
    float total_dist = 0.0f;
    int total_count = 0;

    auto kx = 2.0f / (w - 1);
    auto ky = 2.0f / (h - 1);

    // Compute distance between neighboring pixels in world space,
    // in order to get a good estimate for the initial photon size.
#ifdef USE_STD_THREAD
    parallel_for(0, h / 8 + (h % 8 > 0 ? 1 : 0), [&](size_t ybin) {
        size_t y = ybin * 8;
        float d = 0; int c = 0;
        for (size_t x = 0; x < w; x += 8) {
            Ray rays[4]; Hit hits[4];
            for (int i = 0; i < 4; i++) {
                rays[i] = scene.camera->gen_ray(
                    (x + (i % 2 ? 4 : 0)) * kx - 1.0f,
                    1.0f - (y + (i / 2 ? 4 : 0)) * ky);
                hits[i] = scene.intersect(rays[i]);
            }
            auto eval_distance = [&] (int i, int j) {
                if (hits[i].tri >= 0 && hits[i].tri == hits[j].tri) {
                    d += length((rays[i].org + hits[i].t * rays[i].dir) -
                                (rays[j].org + hits[j].t * rays[j].dir));
                    c++;
                }
            };
            eval_distance(0, 1);
            eval_distance(2, 3);
            eval_distance(0, 2);
            eval_distance(1, 3);
        }

        std::atomic_ref<float>(total_dist) += d;
        std::atomic_ref<int>(total_count) += c;
    });
#else
    #pragma omp parallel for
    for (size_t y = 0; y < h; y += 8) {
        float d = 0; int c = 0;
        for (size_t x = 0; x < w; x += 8) {
            Ray rays[4]; Hit hits[4];
            for (int i = 0; i < 4; i++) {
                rays[i] = scene.camera->gen_ray(
                    (x + (i % 2 ? 4 : 0)) * kx - 1.0f,
                    1.0f - (y + (i / 2 ? 4 : 0)) * ky);
                hits[i] = scene.intersect(rays[i]);
            }
            auto eval_distance = [&] (int i, int j) {
                if (hits[i].tri >= 0 && hits[i].tri == hits[j].tri) {
                    d += length((rays[i].org + hits[i].t * rays[i].dir) -
                                (rays[j].org + hits[j].t * rays[j].dir));
                    c++;
                }
            };
            eval_distance(0, 1);
            eval_distance(2, 3);
            eval_distance(0, 2);
            eval_distance(1, 3);
        }

        #pragma omp atomic
        total_dist += d;
        #pragma omp atomic
        total_count += c;
    }
#endif

    return total_count > 0 ? total_dist / (4 * total_count) : 1.0f;
}

std::unique_ptr<Renderer> create_ppm_renderer(const Scene& scene, size_t max_path_len) {
    return std::unique_ptr<Renderer>(new PhotonMappingRenderer(scene, max_path_len));
}
