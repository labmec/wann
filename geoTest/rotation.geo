// Counterclockwise rotation about z, measured from the positive x-axis.
DefineConstant[well_angle = {90.0, Name "well_angle"}]; // Degrees
rotation_cos = Cos(well_angle * Pi / 180.);
rotation_sin = Sin(well_angle * Pi / 180.);
// Snap axis-aligned rotations to avoid tiny coordinate offsets.
If (Abs(rotation_cos) < 1e-14)
    rotation_cos = 0;
EndIf
If (Abs(rotation_sin) < 1e-14)
    rotation_sin = 0;
EndIf
