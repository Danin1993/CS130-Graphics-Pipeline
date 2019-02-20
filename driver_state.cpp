#include "driver_state.h"
#include <cstring>

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
    std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
    
    state.image_color = new pixel[width * height];
    set_render_black(state);
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
        triangles = state.num_vertices / 3;
        
        for (int i = 0; i < triangles; i++) {
            fill_data_geo(state, &data_geos, vert_index);
            calc_data_geo_pos(state, &data_geos);
            rasterize_triangle(state, (const data_geometry **)(&data_geos));
        }
        break;

    case render_type::indexed:

        break;

    case render_type::fan:

        break;

    case render_type::strip:

        break;

    default:
        std::cerr << "ERROR: invalid render_type specified." << std::endl;
    }

    delete[] data_geos;
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
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
    
    // k0, k1, and k2 are the coefficients for the calculations of the
    // areas for the barycentric coordinates
    float k0[VERT_PER_TRI];
    float k1[VERT_PER_TRI];
    float k2[VERT_PER_TRI]; 
    float total_area;
    float bary[VERT_PER_TRI];
    
    // Calculate pixel coords of vertices
    for (int iter = 0; iter < VERT_PER_TRI; iter++) {
        calc_pixel_coords(state, (*in)[iter], x[iter], y[iter]);
        // Draw pixel at position (i, j)
        //state.image_color[x[iter] + y[iter] * state.image_width] = 
        //    make_pixel(255, 255, 255);
    }

    // Draw the triangle
    // First barycentric weights
    total_area = .5f * ((x[V_B] * y[V_C] - x[V_C] * y[V_B]) 
                        - (x[V_A] * y[V_C] - x[V_C] * y[V_A])
                        - (x[V_A] * y[V_B] - x[V_B] * y[V_A]));

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

    // Iterate through each pixel and calculate the barycentric weights for
    // each.
    for (int y = 0; y < state.image_height; y++) {
        for (int x = 0; x < state.image_width; x++) {
            for (int vert = 0; vert < VERT_PER_TRI; vert++) {
                // Calculation is not done doing the iterative approach
                // We're multiplying every time to find the barycentric
                bary[vert] = .5f * (k0[vert] + (k1[vert] * x) 
                    + (k2[vert] * y)) / total_area;
            }
    
            if (is_pixel_inside(bary)) {
                // At some point this will need to be changed to get the
                // actual color of the pixel.
                state.image_color[x + y * state.image_width] =
                    make_pixel(255, 255, 255);
            }
        }
    }
    

}

void set_render_black(driver_state& state) {
    int image_len = state.image_width * state.image_height;
    for (unsigned i = 0; i < image_len; i++) {
        state.image_color[i] = make_pixel(0, 0, 0);
    }
}

void fill_data_geo(driver_state& state, data_geometry * data_geos[3], 
    int & vert_index) {
    
    for (int i = 0; i < VERT_PER_TRI; i++) {
        (*data_geos)[i].data = state.vertex_data + vert_index;
        vert_index += state.floats_per_vertex;
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
    // i and j might need to be floats, I'm not really sure
    // I mean it works this way, but will it always work?
    i = (int)(w2 * data_geo.gl_Position[X] + (w2 - .5f));
    j = (int)(h2 * data_geo.gl_Position[Y] + (h2 - .5f));
}

bool is_pixel_inside(float * bary_weights) {
    for (int i = 0; i < VERT_PER_TRI; i++) {
        if (bary_weights[i] < 0) {
            return false;
        }
    }

    return true;
}
