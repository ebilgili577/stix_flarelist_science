import argparse

from . import simulator as sm
from scipy.io import readsav
import os
import numpy as np

OUT_FOLDER = 'data/synthetic'

def get_sim_parameters():
    uv_data = readsav('data/stix_data/uv.sav')
    u = uv_data['u']
    v = uv_data['v']

    phase_corr_data = readsav('data/stix_data/phase_corr.sav')
    phase_corr = phase_corr_data['phase_corr']

    pitch_slit_data = readsav('data/stix_data/pitch_slit.sav')
    sf = pitch_slit_data['sf']
    sr = pitch_slit_data['sr']
    pf = pitch_slit_data['pf']
    pr = pitch_slit_data['pr']

    pixel_area = 0.22 * 0.92
    pixel_phase_factor = 45.

    return u, v, sf, sr, pf, pr, phase_corr, pixel_area, pixel_phase_factor

def sim_data(n_samples: int, x_min, x_max, y_min, y_max, add_noise: bool = True , n_sources: int = 1):
    u, v, sf, sr, pf, pr, phase_corr, pixel_area, pixel_phase_factor = get_sim_parameters()
    X = []
    Y = []

    xc = []
    yc = []
    FWHM_max = []
    FWHM_min = []
    flux = []
    orientation = []
    mapcenter = []

    for i in range(n_samples):
        this_xc, this_yc, this_FWHM_max, this_FWHM_min, this_flux, this_orientation, this_mapcenter, this_counts, this_counts_err = sm.SimulateConfig(
            n_sources, u, v, sf, sr, pf, pr, phase_corr, pixel_area, pixel_phase_factor, add_noise=add_noise, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max
            )

        xc.append(this_xc)
        yc.append(this_yc)
        FWHM_max.append(this_FWHM_max)
        FWHM_min.append(this_FWHM_min)
        flux.append(this_flux)
        orientation.append(this_orientation)
        mapcenter.append(this_mapcenter)

        X.append(this_counts)
        Y.append((float(this_xc), float(this_yc)))

    X = np.array(X)
    Y = np.array(Y)
    xc = np.array(xc)
    yc = np.array(yc)
    FWHM_max = np.array(FWHM_max)
    FWHM_min = np.array(FWHM_min)
    flux = np.array(flux)
    orientation = np.array(orientation)
    mapcenter = np.array(mapcenter)
    print(f"Synthetic data X shape is: {X.shape}")
    print(f"Synthetic data Y shape is: {Y.shape}")



    save_data(X=X,Y=Y,xc=xc,yc=yc,FWHM_max=FWHM_max,FWHM_min=FWHM_min,flux=flux,orientation=orientation,mapcenter=mapcenter, n_samples=n_samples, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max)


def save_data(X, Y, xc, yc, FWHM_max, FWHM_min, flux, orientation, mapcenter, n_samples, x_min, x_max, y_min, y_max):
    os.makedirs(OUT_FOLDER, exist_ok=True)
    np.savez(os.path.join(OUT_FOLDER, f'sim_{n_samples}_x{x_min}-{x_max}_y{y_min}-{y_max}.npz'), X=X, Y=Y, xc=xc, yc=yc, FWHM_max=FWHM_max,
             FWHM_min=FWHM_min, flux=flux, orientation=orientation, mapcenter=mapcenter)




def main():
    parser = argparse.ArgumentParser(description="Generate synthetic STIX dataset.")
    parser.add_argument(
        "--n_samples",
        type=int,
        required=True,
        help="Number of simulated samples",
    )
    parser.add_argument(
        "--x_min",
        type=int,
        required=True,
        help="FOV_x min value",
    )
    parser.add_argument(
        "--x_max",
        type=int,
        required=True,
        help="FOV_x max value",
    )
    # TODO: Update doc
    parser.add_argument(
        "--y_min",
        type=int,
        required=True,
        help="FOV_y min value",
    )
    parser.add_argument(
        "--y_max",
        type=int,
        required=True,
        help="FOV_y max value",
    )
    parser.add_argument(
        "--noise",
        type=int,
        choices=[0, 1],
        default=1,
        required=False,
        help="1 = add noise, 0 = no noise (default: 1)",
    )

    args = parser.parse_args()
    add_noise = bool(args.noise)


    filename = f"sim_{args.n_samples}_x{args.x_min}-{args.x_max}_y{args.y_min}-{args.y_max}.npz"
    out_path = os.path.join(OUT_FOLDER, filename)

    if os.path.exists(out_path):
        print(f"Dataset already exists at {out_path}, skipping generation.")
        return

    print(
        f"Simulating data with n_samples={args.n_samples}, "
        f"fov_x_range={args.x_min}-{args.x_max}, "
        f"fov_y_range={args.y_min}-{args.y_max}",
        f"adding noise: {add_noise}",
    )
    sim_data(
        n_samples=args.n_samples,
        x_min=args.x_min,
        x_max=args.x_max,
        y_min=args.y_min,
        y_max=args.y_max,
        add_noise=add_noise,
    )
    print("Simulation finished and data saved.")


if __name__ == "__main__":
    main()