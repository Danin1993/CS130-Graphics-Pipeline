#include "driver_state.h"
#include <cstring>
#include <algorithm>
#include <climits>
#include <cfloat>
#include <vector>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color=0;
    state.image_depth=0;
    
    state.image_len = width * height;

    state.image_color = new pixel[state.image_len];
    set_render_black(state);

    state.image_depth = new float[state.image_len];
    init_image_depth(state);
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    int triangles;
    int vert_index = 0;
    std::cout<<"TODO: implement rendering."<<std::endl;
    
    data_geometry * data_geos = new data_geometry[VERT_PER_TRI];

    switch (type) {
    case render_type::triangle:
        triangles = state.num_vertices / VERT_PER_TRI;
        
        for (int i = 0; i < triangles; i++) {
            fill_data_geo(state, &data_geos, vert_index);
            calc_data_geo_pos(state, &data_geos);
            rasterize_triangle(state, (const data_geometry **)&data_geos);
            //clip_triangle(state, (const data_geometry **)(&data_geos), 0);
        }
        break;

    case render_type::indexed:
        triangles = state.num_vertices / VERT_PER_TRI;
        for (int i = 0; i < triangles; i++) {
            fill_data_geos_indexed(state, &data_geos, vert_index); 
            calc_data_geo_pos(state, &data_geos);
            rasterize_triangle(state, (const data_geometry **)&data_geos);
        }
        break;

    case render_type::fan:

        break;

    case render_type::strip:

        break;

    default:
        std::cerr << "ERROR: Invalid render_type specified." << std::endl;
    }

    delete[] data_geos;
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    std::vector<data_geometry *> tris;
    int sign = 2 * (face % 2) - 1;
    unsigned axis = face % 3;
    bool inside[VERT_PER_TRI] = {0};

    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    
    add_data_geos(state, tris, in);

    for (unsigned i = 0; i < VERT_PER_TRI; i++) {
        // sign will always be either -1 or +1 depending on face
        if (sign > 0) {
            inside[i] = (*in)[i].gl_Position[axis] <
                (*in)[i].gl_Position[W];
        }
        if (sign < 0) {
            inside[i] = (*in)[i].gl_Position[axis] > 
                -1 * (*in)[i].gl_Position[W];
        }
    }

    // The triangle is completely outside of the plane so we don't need to
    // do anything.
    if (all_outside(inside)) { return; }

    // If at least one vertex is inside, then we need to clip
    if (!all_inside(inside) && !all_outside(inside)) {
        // Two vertices outside: (A, B, C) where 1 is inside and 0 outside
        if (inside[0] && !inside[1] && !inside[2]) { // (1, 0, 0)
            create_triangle_2_out(tris, axis, sign, V_A, V_B, V_C, state);
        } else if (!inside[0] && inside[1] && !inside[2]) { // (0, 1, 0)
            create_triangle_2_out(tris, axis, sign, V_B, V_C, V_A, state);
        } else if (!inside[0] && !inside[1] && inside[2]) { // (0, 0, 1)
            create_triangle_2_out(tris, axis, sign, V_C, V_A, V_B, state);
        }
        // One vertex outside: (A, B, C) where 1 is inside and 0 outside
        else if (!inside[0] && inside[1] && inside[2]) { // (0, 1, 1)
            create_triangle_2_in(tris, axis, sign, V_A, V_B, V_C, state);
        } else if (inside[0] && !inside[1] && inside[2]) { // (1, 0, 1)
            create_triangle_2_in(tris, axis, sign, V_B, V_C, V_A, state);
        } else if (inside[0] && inside[1] && !inside[2]) { // (1, 1, 0)
            create_triangle_2_in(tris, axis, sign, V_C, V_A, V_B, state);
        }

        
    }

    // Clip each triangle we've created against the next plane
    for (unsigned i = 0; i < tris.size(); i++) {
        clip_triangle(state,(const data_geometry **)(&(tris[i])),face+1);
    }
    
    clear_data_geos(tris);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    // x and y correspond to the x and y pixel coordinates for each
    // vertex of the triangle.
    int x[VERT_PER_TRI];
    int y[VERT_PER_TRI];

    // z holds the perspective transformed z coordinates for each vertex
    float z[VERT_PER_TRI];
    float depth;

    int min_x, min_y;
    int max_x, max_y;

    
    // k0, k1, and k2 are the coefficients for the calculations of the
    // areas for the barycentric coordinates
    float k0[VERT_PER_TRI];
    float k1[VERT_PER_TRI];
    float k2[VERT_PER_TRI]; 
    float total_area;
    float bary[VERT_PER_TRI];

    // Define and alloc the data_frag for later
    data_fragment frag;
    // Allocate an array to hold the interpolated data
    frag.data = new float[MAX_FLOATS_PER_VERTEX];
    

    // Calculate pixel coords of vertices
    for (int iter = 0; iter < VERT_PER_TRI; iter++) {
        // Conversion to homogenous coords (for x and y) is done in this
        // function. Do not forget to do it for z.
        calc_pixel_coords(state, (*in)[iter], x[iter], y[iter]);
        //std::cout << "DEBUG: (" << x[iter] << ", " << y[iter] << ")\n";
    }

    // Draw the triangle
    // First barycentric weights
    total_area = .5f * ((x[V_B] * y[V_C] - x[V_C] * y[V_B]) 
                        - (x[V_A] * y[V_C] - x[V_C] * y[V_A])
                        + (x[V_A] * y[V_B] - x[V_B] * y[V_A]));

    // These are constant for all pixels so let's calculate them ahead of
    // time.
    k0[V_A] = x[V_B] * y[V_C] - x[V_C] * y[V_B];
    k1[V_A] = y[V_B] - y[V_C];
    k2[V_A] = x[V_C] - x[V_B];

    k0[V_B] = x[V_C] * y[V_A] - x[V_A] * y[V_C];
    k1[V_B] = y[V_C] - y[V_A];
    k2[V_B] = x[V_A] - x[V_C];

    k0[V_C] = x[V_A] * y[V_B] - x[V_B] * y[V_A];
    k1[V_C] = y[V_A] - y[V_B];
    k2[V_C] = x[V_B] - x[V_A];

    // Calculate the min and max coordinates to iterate through
    calc_min_coord(state, x, y, min_x, min_y);
    calc_max_coord(state, x, y, max_x, max_y);

    /*
    std::cout << "DEBUG: min: (" << min_x << ", " << min_y << ")\n";
    std::cout << "DEBUG: max: (" << max_x << ", " << max_y << ")";
    */

    // Iterate through each pixel and calculate the barycentric weights for
    // each.
    for (int y = min_y + 1; y < max_y + 1; y++) {
        for (int x = min_x + 1; x < max_x + 1; x++) {
            for (int vert = 0; vert < VERT_PER_TRI; vert++) {
                // Calculation is not done doing the iterative approach
                // We're multiplying every time to find the barycentric
                bary[vert] = .5f * (k0[vert] + (k1[vert] * x) 
                    + (k2[vert] * y)) / total_area;
            }


            calc_z_coords(in, z);
            depth = calc_depth_at(z, bary);

            // Only draw if the pixel is inside the triangle and it is the
            // closest triangle to the camera
            if (is_pixel_inside(bary) && 
                depth < state.image_depth[x + y * state.image_width]) {

                state.image_color[x + y * state.image_width] =
                    get_pixel_color(state, frag, in, bary);
                state.image_depth[x + y * state.image_width] = depth;
            }
        }
    }

    // Don't forget to delete our allocated memory!  
    delete[] frag.data;

}


/**************************************************************************/
/* Initialization */
/**************************************************************************/

void set_render_black(driver_state& state) {
    for (int i = 0; i < state.image_len; i++) {
        state.image_color[i] = make_pixel(0, 0, 0);
    }
}

void init_image_depth(driver_state& state) {
    for (int i = 0; i < state.image_len; i++) {
        state.image_depth[i] = FLT_MAX;
    }
}


/**************************************************************************/
/* Rasterize Triangle Helpers */
/**************************************************************************/

void fill_data_geo(driver_state& state, data_geometry * data_geos[3], 
    int & vert_index) {
    
    for (int i = 0; i < VERT_PER_TRI; i++) {
        (*data_geos)[i].data = state.vertex_data + vert_index;
        vert_index += state.floats_per_vertex;
    }
    
}

void fill_data_geos_indexed(driver_state& state,
    data_geometry * data_geos[3], int & vert_index) {
    
    for (int i = 0; i < VERT_PER_TRI; i++) {
        (*data_geos)[i].data = state.vertex_data 
            + state.index_data[vert_index] * state.floats_per_vertex;
        vert_index += VERT_PER_TRI;
    }
}

void calc_data_geo_pos(driver_state& state, data_geometry * data_geos[3]) {
    data_vertex data_vert;
    for (int i = 0; i < VERT_PER_TRI; i++) {
        data_vert.data = (*data_geos)[i].data;
        state.vertex_shader(data_vert, (*data_geos)[i], state.uniform_data);
    }
}

void calc_pixel_coords(driver_state& state, const data_geometry& data_geo, 
    int& i, int& j) {
    
    static const float w2 = state.image_width / 2.0f;
    static const float h2 = state.image_height / 2.0f;
    
    // The conversion to homogeneous coords happens here.
    i = (int)(w2 * data_geo.gl_Position[X] / data_geo.gl_Position[W]
        + (w2 - .5f));
    j = (int)(h2 * data_geo.gl_Position[Y] / data_geo.gl_Position[W]
        + (h2 - .5f));
}

void calc_min_coord(const driver_state& state, int * x, int * y, int& min_x,
     int& min_y) {
    
    // The maximum pixel coord we can have is (width - 1, height - 1), so
    // set the starting values to those
    min_x = state.image_width - 1;
    min_y = state.image_height - 1;

    for (int i = 0; i < VERT_PER_TRI; i++) {
        min_x = std::min(min_x, x[i]);
        min_y = std::min(min_y, y[i]);
    }

    // The minimum pixel coord we can have is (0, 0) so if either min_x or 
    // min_y are negative we set them to 0.
    min_x = std::max(min_x, 0);
    min_y = std::max(min_y, 0);
}

void calc_max_coord(const driver_state& state, int * x, int * y, int& max_x,
    int& max_y) {
    
    // The minimum pixel coord we can have is (0, 0), so set the starting
    // value to 0, 0; 
    max_x = 0;
    max_y = 0;

    for (int i = 0; i < VERT_PER_TRI; i++) {
        max_x = std::max(max_x, x[i]);
        max_y = std::max(max_y, y[i]);
    }

    // The maximum pixel coord we can have is (width, height) so if either 
    // min_x or max_y are greater then we set them accordingly.
    max_x = std::min(max_x, state.image_width - 1);
    max_y = std::min(max_y, state.image_height - 1);
}

bool is_pixel_inside(float * bary_weights) {
    for (int i = 0; i < VERT_PER_TRI; i++) {
        if (bary_weights[i] < 0) {
            return false;
        }
    }

    return true;
}


/**************************************************************************/
/* Fragment Shader */
/**************************************************************************/

pixel get_pixel_color(driver_state& state, data_fragment& frag,
    const data_geometry * data_geos[3], float * screen_bary) {
    
    float world_bary[VERT_PER_TRI];
    data_output out;

    // For each float in the vertex we have to interpolate data depending
    // on the interp_rule associated with it.
    for (int i = 0; i < state.floats_per_vertex; i++) {
        switch (state.interp_rules[i]) {
        
        // If the interpolation rule is flat then set all data floats equal
        // to the data of the first vertex
        case interp_type::flat:
            frag.data[i] = (*data_geos)[V_A].data[i];
            break;

        // If the interpolation rule is smooth then we want perspective
        // correct interpolation
        case interp_type::smooth:
            // Convert the bary centric coordinates we calculated previously
            // in screen space back to world space
            convert_from_screen(screen_bary, world_bary, data_geos);
            frag.data[i] = interpolate_fragment_at(i, data_geos, 
                world_bary);
            break;

        // If the interpolation rule is noperspective then we just want
        // interpolation based on our screen space barycentric coordinates
        case interp_type::noperspective:
            frag.data[i] = interpolate_fragment_at(i, data_geos,
                screen_bary);
            break;

        default:
            std::cerr << "ERROR: Invalid interp_type specified.\n";
            break;
        }
    }

    // Call our fragment shader with the data we just interpolated
    state.fragment_shader(frag, out, state.uniform_data);

    // Multiply the output by C_MAX (255) because output_color returns a
    // value [0, 1]
    return make_pixel(out.output_color[C_R] * C_MAX, out.output_color[C_G] 
        * C_MAX, out.output_color[C_B] * C_MAX);
}


float interpolate_fragment_at(unsigned index,
    const data_geometry * data_geos[3], float * bary) {
    float ret = 0;

    for (unsigned i = 0; i < VERT_PER_TRI; i++) {
        ret += bary[i] * (*data_geos)[i].data[index];
    }

    return ret;
}

void convert_from_screen(float * screen_bary, float * world_bary, 
    const data_geometry * data_geos[3]) {
    float k = 0;

    for (unsigned i = 0; i < VERT_PER_TRI; i++) {
        k += screen_bary[i] / (*data_geos)[i].gl_Position[W];
    }

    for (unsigned i = 0; i < VERT_PER_TRI; i++) {
        world_bary[i] = screen_bary[i] / ((*data_geos)[i].gl_Position[W]
            * k);
    }
}


/**************************************************************************/
/* Z-Buffer */
/**************************************************************************/

void calc_z_coords(const data_geometry * data_geos[3], float * z) {
    
    for (unsigned i = 0; i < VERT_PER_TRI; i++) {
        z[i] = (*data_geos)[i].gl_Position[Z] 
            / (*data_geos)[i].gl_Position[W];
    }
}

float calc_depth_at(float * z, float * bary) {
    float ret = 0;

    for (unsigned i = 0; i < VERT_PER_TRI; i++) {
        ret += z[i] * bary[i];
    }

    return ret;
}


/**************************************************************************/
/* Clipping */
/**************************************************************************/
void clear_data_geos(std::vector<data_geometry *>& tris) {
    for (unsigned i = 0; i < tris.size(); i++) {
        remove_data_geo(i, tris);
    }

    tris.clear();
}

void remove_data_geo(unsigned index, std::vector<data_geometry *>& tris) {
    for (unsigned i = 0; i < VERT_PER_TRI; i++) {
        delete[] tris[index][i].data;
    }
    delete[] tris[index];
    tris.erase(tris.begin() + index);
}

void add_data_geos(const driver_state& state, 
    std::vector<data_geometry *>& tris, 
    const data_geometry * data_geos[3]) {

    tris.push_back(new data_geometry[VERT_PER_TRI]);
    for (unsigned i = 0; i < VERT_PER_TRI; i++) {
        tris[tris.size() - 1][i].data = new float[MAX_FLOATS_PER_VERTEX];
        for (unsigned j = 0; j < DATA_PER_COORD; j++) {
            tris[tris.size() - 1][i].gl_Position[j] =
                (*data_geos)[i].gl_Position[j];
        }

        for (int j = 0; j < state.floats_per_vertex; j++) {
            tris[tris.size() - 1][i].data[j] = tris[0][i].data[j];
        }
    }
}

void add_data_geos(std::vector<data_geometry *>& tris, vec4 a, vec4 b,
    vec4 c) {

    data_geometry * geos = new data_geometry[VERT_PER_TRI];

    for (unsigned i = 0; i < VERT_PER_TRI; i++) {
        geos[i].data = new float[MAX_FLOATS_PER_VERTEX];
    }

    geos[V_A].gl_Position = a;
    geos[V_B].gl_Position = b;
    geos[V_C].gl_Position = c;

    tris.push_back(geos);
}

void copy_data_geos_data(const driver_state& state,
    const data_geometry& from, data_geometry& to) {

    for (int i = 0; i < state.floats_per_vertex; i++) {
        to.data[i] = from.data[i];
    }
}

bool all_inside(bool * inside) {
    for (unsigned i = 0; i < VERT_PER_TRI; i++) {
        if (!inside[i]) {
            return false;
        }
    }

    return true;
}

bool all_outside(bool * inside) {
    for (unsigned i = 0; i < VERT_PER_TRI; i++) {
        if (inside[i]) {
            return false;
        }
    }

    return true;
}

void create_triangle_2_out(std::vector<data_geometry *>& tris,
    unsigned axis, int sign, unsigned in_index, unsigned out0_index,
    unsigned out1_index, const driver_state& state) {
    
    data_geometry a = tris[0][in_index];
    data_geometry b = tris[0][out0_index];
    data_geometry c = tris[0][out1_index];

    float b_numer = sign * b.gl_Position[W] - b.gl_Position[axis];
    float ab_weight = b_numer / ((a.gl_Position[axis]
        - sign * a.gl_Position[W]) + b_numer);
    float c_numer = sign * c.gl_Position[W] - c.gl_Position[axis];
    float ac_weight = c_numer / ((a.gl_Position[axis] 
        - sign * a.gl_Position[W]) + c_numer);

    float ab_nop_weight;
    float ac_nop_weight;

    add_data_geos(tris, a.gl_Position, ab_weight * a.gl_Position 
        + (1 - ab_weight) * b.gl_Position, ac_weight * a.gl_Position 
        + (1 - ac_weight) * c.gl_Position);


    // All interp_types copy the inside index to the a vertex
    copy_data_geos_data(state, tris[0][in_index], tris[1][V_A]);
    
    for (int i = 0; i < state.floats_per_vertex; i++) {
        switch (state.interp_rules[i]) {
        case interp_type::flat:
            for (unsigned j = 1; j < VERT_PER_TRI; j++) {
                tris[1][j].data[i] = tris[0][in_index].data[i];
            }
            break;
    
        case interp_type::smooth:
            tris[1][V_B].data[i] = interpolate_data(ab_weight,
                tris[0][in_index].data[i], tris[0][out0_index].data[i]);
                
            tris[1][V_C].data[i] = interpolate_data(ac_weight,
                tris[0][in_index].data[i], tris[0][out1_index].data[i]);
            break;

        case interp_type::noperspective:
            ab_nop_weight = ab_weight * tris[0][in_index].gl_Position[W]
                * calc_noperspective_weight(ab_weight, 
                    tris[0][in_index].gl_Position[W], 
                    tris[0][out0_index].gl_Position[W]);

            ac_nop_weight = ac_weight * tris[0][in_index].gl_Position[W]
                * calc_noperspective_weight(ac_weight,
                    tris[0][in_index].gl_Position[W],
                    tris[0][out1_index].gl_Position[W]);

            tris[1][V_B].data[i] = interpolate_data(ab_nop_weight,
                tris[0][in_index].data[i], tris[0][out0_index].data[i]);
                
            tris[1][V_C].data[i] = interpolate_data(ac_nop_weight,
                tris[0][in_index].data[i], tris[0][out1_index].data[i]);
            break;

        default:
            std::cerr << "ERROR: Invalid interp_type specified.\n";
            break;
        }
    }
    remove_data_geo(0, tris);
}

void create_triangle_2_in(std::vector<data_geometry *>& tris,
    unsigned axis, int sign, unsigned out_index, unsigned in0_index,
    unsigned in1_index, const driver_state& state) {

    
    data_geometry a = tris[0][out_index];
    data_geometry b = tris[0][in0_index];
    data_geometry c = tris[0][in1_index];

    float b_numer = sign * b.gl_Position[W] - b.gl_Position[axis];
    float ab_weight = b_numer / ((a.gl_Position[axis]
        - sign * a.gl_Position[W]) + b_numer);
    float c_numer = sign * c.gl_Position[W] - c.gl_Position[axis];
    float ac_weight = c_numer / ((a.gl_Position[axis] 
        - sign * a.gl_Position[W]) + c_numer);

    float ab_nop_weight;
    float ac_nop_weight;

/*
    add_data_geos(tris, b.gl_Position, c.gl_Position, 
        ac_weight * a.gl_Position + (1 - ac_weight) * c.gl_Position);
    add_data_geos(tris, ac_weight * a.gl_Position 
        + (1 - ac_weight) * c.gl_Position, ab_weight * a.gl_Position 
        + (1 - ab_weight) * b.gl_Position, b.gl_Position);
*/

    add_data_geos(state, tris, (const data_geometry **)(&tris[0]));
    tris[tris.size() - 1][V_A].gl_Position = 
        ac_weight * a.gl_Position + (1 - ac_weight) * c.gl_Position;

    add_data_geos(state, tris, (const data_geometry **)(&tris[0]));
    tris[tris.size() - 1][V_A].gl_Position =
        ab_weight * a.gl_Position + (1 - ab_weight) * b.gl_Position;
    tris[tris.size() - 1][V_C].gl_Position = tris[tris.size() - 2][V_A].gl_Position;
        
        
    for (int i = 0; i < state.floats_per_vertex; i++) {
        switch (state.interp_rules[i]) {
        case interp_type::flat:
            for (unsigned j = 1; j < VERT_PER_TRI; j++) {
                tris[1][j].data[i] = tris[0][out_index].data[i];
                tris[2][j].data[i] = tris[0][out_index].data[i];
            }
            break;
    
        case interp_type::smooth:
            tris[1][V_A].data[i] = interpolate_data(ac_weight,
                tris[0][out_index].data[i], tris[0][in1_index].data[i]);
            
            copy_data_geos_data(state, tris[1][V_A], tris[2][V_C]);
             
            tris[2][V_A].data[i] = interpolate_data(ab_weight,
                tris[0][out_index].data[i], tris[0][in0_index].data[i]);
            break;

        case interp_type::noperspective:
            ab_nop_weight = ab_weight * tris[0][out_index].gl_Position[W]
                * calc_noperspective_weight(ab_weight, 
                    tris[0][out_index].gl_Position[W], 
                    tris[0][in0_index].gl_Position[W]);

            ac_nop_weight = ac_weight * tris[0][out_index].gl_Position[W]
                * calc_noperspective_weight(ac_weight,
                    tris[0][out_index].gl_Position[W],
                    tris[0][in1_index].gl_Position[W]);

            tris[1][V_A].data[i] = interpolate_data(ac_nop_weight,
                tris[0][out_index].data[i], tris[0][in1_index].data[i]);
                
            copy_data_geos_data(state, tris[1][V_A], tris[2][V_C]); 
            tris[2][V_A].data[i] = interpolate_data(ab_nop_weight,
                tris[0][out_index].data[i], tris[0][in0_index].data[i]);
            break;

        default:
            std::cerr << "ERROR: Invalid interp_type specified.\n";
            break;
        }
    }
    remove_data_geo(0, tris);
}

float interpolate_data(float weight, float data0, float data1) {
    return weight * data0 + (1 - weight) * data1;
}

float calc_noperspective_weight(float weight, float a_w, float p_w) {
    return 1.0f / (weight * a_w + (1 - weight) * p_w);
}
