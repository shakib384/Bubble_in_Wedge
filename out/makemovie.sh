#!/bin/bash

ffmpeg -y -framerate 40 -pattern_type glob -i 'vof-omega_y_view1_*.png' -c:v libx264 -crf 18 -pix_fmt yuv420p vof-omega_y_view1.mp4
ffmpeg -y -framerate 40 -pattern_type glob -i 'vof-p_view1_*.png' -c:v libx264 -crf 18 -pix_fmt yuv420p vof-p_view1.mp4
ffmpeg -y -framerate 40 -pattern_type glob -i 'vof-p_view2_*.png' -c:v libx264 -crf 18 -pix_fmt yuv420p vof-p_view2.mp4
ffmpeg -y -framerate 40 -pattern_type glob -i 'vof-ux_view1_*.png' -c:v libx264 -crf 18 -pix_fmt yuv420p vof-ux_view1.mp4
ffmpeg -y -framerate 40 -pattern_type glob -i 'vof-uz_view1_*.png' -c:v libx264 -crf 18 -pix_fmt yuv420p vof-uz_view1.mp4
ffmpeg -y -framerate 40 -pattern_type glob -i 'vof-uy_view4_*.png' -c:v libx264 -crf 18 -pix_fmt yuv420p vof-uy_view4.mp4
ffmpeg -y -framerate 40 -pattern_type glob -i 'vof-uz_view4_*.png' -c:v libx264 -crf 18 -pix_fmt yuv420p vof-uz_view4.mp4
ffmpeg -y -framerate 40 -pattern_type glob -i 'vof-p_view4_*.png' -c:v libx264 -crf 18 -pix_fmt yuv420p vof-p_view4.mp4
