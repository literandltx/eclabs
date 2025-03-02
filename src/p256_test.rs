use num_bigint::BigInt;
use num_traits::{Num, One};

#[cfg(test)]
mod correctness {
    use super::*;
    use crate::{helpers, Actor, AffinePoint, Curve, ProjectivePoint};

    #[test]
    fn test_map_to_affine_point_trivial() {
        let p: BigInt = BigInt::from_str_radix("26", 10).unwrap();
        let a: BigInt = BigInt::from_str_radix("1", 10).unwrap();
        let b: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let curve: Curve = Curve::new(p, a, b).unwrap();
        let x: BigInt = BigInt::from_str_radix("2", 16).unwrap();
        let y: BigInt = BigInt::from_str_radix("3", 16).unwrap();

        // (2, 3, 1)
        let projective_point: ProjectivePoint =
            ProjectivePoint::new(x.clone(), y.clone(), BigInt::one());

        // (2, 3)
        let affine_point: AffinePoint = projective_point.to_affine(&curve).unwrap();

        assert_eq!(&affine_point.x, &x);
        assert_eq!(&affine_point.y, &y);
    }

    #[test]
    fn test_no_affine_form_of_neutral() {
        let p: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
            16,
        )
        .unwrap();
        let a: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000fffffffffffffffffffffffc",
            16,
        )
        .unwrap();
        let b: BigInt = BigInt::from_str_radix(
            "5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
            16,
        )
        .unwrap();

        let curve: Curve = Curve::new(p, a, b).unwrap();

        let x1: BigInt = BigInt::from_str_radix("0", 10).unwrap();
        let y1: BigInt = BigInt::from_str_radix("1", 10).unwrap();
        let z1: BigInt = BigInt::from_str_radix("0", 10).unwrap();

        let point1: ProjectivePoint = ProjectivePoint::new(x1.clone(), y1.clone(), z1.clone());

        assert!(point1.to_affine(&curve).is_none());
    }

    #[test]
    fn test_map_to_affine_point_inverse_exist() {
        let p: BigInt = BigInt::from_str_radix("26", 10).unwrap();
        let a: BigInt = BigInt::from_str_radix("1", 10).unwrap();
        let b: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let curve: Curve = Curve::new(p, a, b).unwrap();
        let x: BigInt = BigInt::from_str_radix("2", 10).unwrap();
        let y: BigInt = BigInt::from_str_radix("3", 10).unwrap();
        let z: BigInt = BigInt::from_str_radix("11", 10).unwrap();

        // (2, 3, 1)
        let projective_point: ProjectivePoint =
            ProjectivePoint::new(x.clone(), y.clone(), z.clone());

        // (12, 5)
        let affine_point: AffinePoint = projective_point.to_affine(&curve).unwrap();

        assert_eq!(affine_point.x, BigInt::from(12));
        assert_eq!(affine_point.y, BigInt::from(5));
    }

    #[test]
    fn test_map_to_affine_point_inverse_non_exist() {
        let p: BigInt = BigInt::from_str_radix("26", 10).unwrap();
        let a: BigInt = BigInt::from_str_radix("1", 10).unwrap();
        let b: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let curve: Curve = Curve::new(p, a, b).unwrap();
        let x: BigInt = BigInt::from_str_radix("2", 10).unwrap();
        let y: BigInt = BigInt::from_str_radix("3", 10).unwrap();
        let z: BigInt = BigInt::from_str_radix("10", 10).unwrap();

        // (2, 3, 1)
        let projective_point: ProjectivePoint =
            ProjectivePoint::new(x.clone(), y.clone(), z.clone());

        // (12, 5)
        let affine_point: Option<AffinePoint> = projective_point.to_affine(&curve);
        assert!(&affine_point.is_none());
    }

    #[test]
    fn test_double1() {
        let p: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
            16,
        )
        .unwrap();
        let a: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000fffffffffffffffffffffffc",
            16,
        )
        .unwrap();
        let b: BigInt = BigInt::from_str_radix(
            "5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
            16,
        )
        .unwrap();

        let curve: Curve = Curve::new(p, a, b).unwrap();
        let x: BigInt = BigInt::from_str_radix(
            "96f43acb7e4ad862d547681483abfe1c6a21eace101ceca1e0dfa97beb8d99ec",
            16,
        )
        .unwrap();
        let y: BigInt = BigInt::from_str_radix(
            "69e14d9fe6ce5b8cb04b2fb0ffa525526c2b7114d891b4aa2f53dfe8de3d6dad",
            16,
        )
        .unwrap();
        let z: BigInt = BigInt::from_str_radix("1", 16).unwrap();

        let double_x: BigInt = BigInt::from_str_radix(
            "23256713097452873534819684224181198488753197392778987539588939509885686328462",
            10,
        )
        .unwrap();
        let double_y: BigInt = BigInt::from_str_radix(
            "23560293084035690959730279798706588908809082944968261336868665854561491207411",
            10,
        )
        .unwrap();
        let double_z: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let projective_point: ProjectivePoint =
            ProjectivePoint::new(x.clone(), y.clone(), z.clone());

        let double_projective_point: ProjectivePoint = curve.double(&projective_point);
        let double_affine_point = double_projective_point.to_affine(&curve).unwrap();

        let result: ProjectivePoint = double_affine_point.to_projective();

        let tmp = curve.verify_projective_point(&result);
        // println!("{:?}", tmp);

        assert_eq!(result.z, double_z);
        assert_eq!(result.x, double_x);
        assert_eq!(result.y, double_y);
    }

    #[test]
    fn test_add_affine_point_neutral() {
        let p: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
            16,
        )
        .unwrap();
        let a: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000fffffffffffffffffffffffc",
            16,
        )
        .unwrap();
        let b: BigInt = BigInt::from_str_radix(
            "5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
            16,
        )
        .unwrap();

        let curve: Curve = Curve::new(p, a, b).unwrap();

        let x1: BigInt = BigInt::from_str_radix(
            "101762592313467086577367687451288259444322615470288307860676415566273649512588",
            10,
        )
        .unwrap();
        let y1: BigInt = BigInt::from_str_radix(
            "102965935123701046296389518845320850966553667754075485002495435140169891861431",
            10,
        )
        .unwrap();
        let z1: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let point1: ProjectivePoint = ProjectivePoint::new(x1.clone(), y1.clone(), z1.clone());
        let identity: ProjectivePoint = ProjectivePoint::neutral();

        let add_projective_point: ProjectivePoint = curve.add(&point1, &identity);
        let add_affine_point: AffinePoint = add_projective_point.to_affine(&curve).unwrap();
        let result: ProjectivePoint = add_affine_point.to_projective();

        assert_eq!(result.x, x1);
        assert_eq!(result.y, y1);
        assert_eq!(result.z, z1);
    }

    #[test]
    fn test_add_two_affine_points() {
        let p: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
            16,
        )
        .unwrap();
        let a: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000fffffffffffffffffffffffc",
            16,
        )
        .unwrap();
        let b: BigInt = BigInt::from_str_radix(
            "5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
            16,
        )
        .unwrap();

        let curve: Curve = Curve::new(p, a, b).unwrap();

        let x1: BigInt = BigInt::from_str_radix(
            "68278443758753309620248260190700068206105694628838492838163118930873086220780",
            10,
        )
        .unwrap();
        let y1: BigInt = BigInt::from_str_radix(
            "47890925436225047886667813019233725009834729938480004609112948108082127138221",
            10,
        )
        .unwrap();
        let z1: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let x2: BigInt = BigInt::from_str_radix(
            "36041773901858648762789006722537555576806908678254251775422211831047989597049",
            10,
        )
        .unwrap();
        let y2: BigInt = BigInt::from_str_radix(
            "104374463647532952902541679499731110179928673815908965998191930674877316391086",
            10,
        )
        .unwrap();
        let z2: BigInt = BigInt::from_str_radix("1", 16).unwrap();

        let xs: BigInt = BigInt::from_str_radix(
            "60925061547489862940527908037003633368921254513905133864735524490867813970614",
            10,
        )
        .unwrap();
        let ys: BigInt = BigInt::from_str_radix(
            "4224995115825146765508233942604852388022531386580408859355570965592456321292",
            10,
        )
        .unwrap();
        let zs: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let point1: ProjectivePoint = ProjectivePoint::new(x1.clone(), y1.clone(), z1.clone());
        let point2: ProjectivePoint = ProjectivePoint::new(x2.clone(), y2.clone(), z2.clone());

        let add_projective_point: ProjectivePoint = curve.add(&point1, &point2);

        let add_affine_point: AffinePoint = add_projective_point.to_affine(&curve).unwrap();

        let result: ProjectivePoint = add_affine_point.to_projective();

        assert_eq!(result.x, xs);
        assert_eq!(result.y, ys);
        assert_eq!(result.z, zs);
    }

    #[test]
    fn test_multiply_point_on_scalar() {
        let p: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
            16,
        )
        .unwrap();
        let a: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000fffffffffffffffffffffffc",
            16,
        )
        .unwrap();
        let b: BigInt = BigInt::from_str_radix(
            "5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
            16,
        )
        .unwrap();

        let curve: Curve = Curve::new(p, a, b).unwrap();

        let x: BigInt = BigInt::from_str_radix(
            "36468809527719095167323546128622644075555946232974359005635481469739685065250",
            10,
        )
        .unwrap();
        let y: BigInt = BigInt::from_str_radix(
            "55086920727300960012414320103293286937832915861215179864890281344876589042713",
            10,
        )
        .unwrap();
        let z: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let scalar: BigInt = BigInt::from_str_radix(
            "50313277712876843208627822877019546215941224990452903113298304423968512062362",
            10,
        )
        .unwrap();

        let x_res: BigInt = BigInt::from_str_radix(
            "105471300999999635061948442487927318986215238139837487327895353516078990089204",
            10,
        )
        .unwrap();
        let y_res: BigInt = BigInt::from_str_radix(
            "102129873108155434225483148806386549747776446537842754889946033852013914809773",
            10,
        )
        .unwrap();
        let z_res: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let point: ProjectivePoint = ProjectivePoint::new(x.clone(), y.clone(), z.clone());

        let mul_projective_point: ProjectivePoint = curve.scalar_mul(&scalar, &point);

        let mul_affine_point: AffinePoint = mul_projective_point.to_affine(&curve).unwrap();

        let result: ProjectivePoint = mul_affine_point.to_projective();

        assert_eq!(result.x, x_res);
        assert_eq!(result.y, y_res);
        assert_eq!(result.z, z_res);
    }

    #[test]
    fn test_multiply_point_on_scalar_montgomery() {
        let p: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
            16,
        )
        .unwrap();
        let a: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000fffffffffffffffffffffffc",
            16,
        )
        .unwrap();
        let b: BigInt = BigInt::from_str_radix(
            "5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
            16,
        )
        .unwrap();

        let curve: Curve = Curve::new(p, a, b).unwrap();

        let x: BigInt = BigInt::from_str_radix(
            "36468809527719095167323546128622644075555946232974359005635481469739685065250",
            10,
        )
        .unwrap();
        let y: BigInt = BigInt::from_str_radix(
            "55086920727300960012414320103293286937832915861215179864890281344876589042713",
            10,
        )
        .unwrap();
        let z: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let scalar: BigInt = BigInt::from_str_radix(
            "50313277712876843208627822877019546215941224990452903113298304423968512062362",
            10,
        )
        .unwrap();

        let x_res: BigInt = BigInt::from_str_radix(
            "105471300999999635061948442487927318986215238139837487327895353516078990089204",
            10,
        )
        .unwrap();
        let y_res: BigInt = BigInt::from_str_radix(
            "102129873108155434225483148806386549747776446537842754889946033852013914809773",
            10,
        )
        .unwrap();
        let z_res: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let point: ProjectivePoint = ProjectivePoint::new(x.clone(), y.clone(), z.clone());

        let mul_projective_point: ProjectivePoint = curve.scalar_mul_montgomery(&scalar, &point);

        let mul_affine_point: AffinePoint = mul_projective_point.to_affine(&curve).unwrap();

        let result: ProjectivePoint = mul_affine_point.to_projective();

        assert_eq!(result.x, x_res);
        assert_eq!(result.y, y_res);
        assert_eq!(result.z, z_res);
    }

    #[test]
    fn test_scalar_multiplication_using_order() {
        let p: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
            16,
        )
        .unwrap();
        let a: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000fffffffffffffffffffffffc",
            16,
        )
        .unwrap();
        let b: BigInt = BigInt::from_str_radix(
            "5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
            16,
        )
        .unwrap();

        let curve: Curve = Curve::new(p, a, b).unwrap();

        let x: BigInt = BigInt::from_str_radix(
            "88830848574649868340568578826653665181106303545668410351759331370193752307474",
            10,
        )
        .unwrap();
        let y: BigInt = BigInt::from_str_radix(
            "22058354116588620035967384506492701987686288267863232015425921192621894353676",
            10,
        )
        .unwrap();
        let z: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let x_res: BigInt = BigInt::from_str_radix("0", 10).unwrap();
        let y_res: BigInt = BigInt::from_str_radix("1", 10).unwrap();
        let z_res: BigInt = BigInt::from_str_radix("0", 10).unwrap();

        let point: ProjectivePoint = ProjectivePoint::new(x.clone(), y.clone(), z.clone());

        let result: ProjectivePoint = curve.scalar_mul(&curve.get_order(), &point);

        assert_eq!(result.x, x_res);
        assert_eq!(result.y, y_res);
        assert_eq!(result.z, z_res);
    }

    #[test]
    fn test_scalar_multiplication_using_order_montgomery() {
        let p: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
            16,
        )
        .unwrap();
        let a: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000fffffffffffffffffffffffc",
            16,
        )
        .unwrap();
        let b: BigInt = BigInt::from_str_radix(
            "5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
            16,
        )
        .unwrap();

        let curve: Curve = Curve::new(p, a, b).unwrap();

        let x: BigInt = BigInt::from_str_radix(
            "88830848574649868340568578826653665181106303545668410351759331370193752307474",
            10,
        )
        .unwrap();
        let y: BigInt = BigInt::from_str_radix(
            "22058354116588620035967384506492701987686288267863232015425921192621894353676",
            10,
        )
        .unwrap();
        let z: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let x_res: BigInt = BigInt::from_str_radix("0", 10).unwrap();
        let y_res: BigInt = BigInt::from_str_radix("1", 10).unwrap();
        let z_res: BigInt = BigInt::from_str_radix("0", 10).unwrap();

        let point: ProjectivePoint = ProjectivePoint::new(x.clone(), y.clone(), z.clone());

        let result: ProjectivePoint = curve.scalar_mul_montgomery(&curve.get_order(), &point);

        assert_eq!(result.x, x_res);
        assert_eq!(result.y, y_res);
        assert_eq!(result.z, z_res);

        assert!(curve.verify_projective_point(&result));
    }

    #[test]
    fn test_complex_scalar_multiplication_and_sum() {
        let p: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
            16,
        )
        .unwrap();
        let a: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000fffffffffffffffffffffffc",
            16,
        )
        .unwrap();
        let b: BigInt = BigInt::from_str_radix(
            "5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
            16,
        )
        .unwrap();

        let curve: Curve = Curve::new(p, a, b).unwrap();

        let x1: BigInt = BigInt::from_str_radix(
            "29424756536908666275014196764752012937549470486181629932252170285403408049517",
            10,
        )
        .unwrap();
        let y1: BigInt = BigInt::from_str_radix(
            "99699357898675284307950095025727949597877550412833207843633386200194164890561",
            10,
        )
        .unwrap();
        let z1: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let x2: BigInt = BigInt::from_str_radix(
            "95387701464026858343496773958108763791065521288270402154042613065295218611790",
            10,
        )
        .unwrap();
        let y2: BigInt = BigInt::from_str_radix(
            "56197753972877288816938393898610332375031967031226761314470230286071453671148",
            10,
        )
        .unwrap();
        let z2: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let scalar1: BigInt = BigInt::from_str_radix(
            "115166654073940935621181939619999211512712171442155919412146477633437056010824",
            10,
        )
        .unwrap();
        let scalar2: BigInt = BigInt::from_str_radix(
            "18076504642849865357726580190900882117966238578277710228703901228455155726225",
            10,
        )
        .unwrap();

        let x_res: BigInt = BigInt::from_str_radix(
            "37352616003680990833471240280964019197164191581612728497391994362995806189072",
            10,
        )
        .unwrap();
        let y_res: BigInt = BigInt::from_str_radix(
            "51702340847344877490410640050738219062093639452275966547398163436084404331575",
            10,
        )
        .unwrap();
        let z_res: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let point1: ProjectivePoint = ProjectivePoint::new(x1.clone(), y1.clone(), z1.clone());
        let point2: ProjectivePoint = ProjectivePoint::new(x2.clone(), y2.clone(), z2.clone());

        let computed_result_projective_point: ProjectivePoint = curve.add(
            &curve.scalar_mul(&scalar1, &point1),
            &curve.scalar_mul(&scalar2, &point2),
        );

        let computed_result_affine_point: AffinePoint =
            computed_result_projective_point.to_affine(&curve).unwrap();
        let result: ProjectivePoint = computed_result_affine_point.to_projective();

        assert_eq!(result.x, x_res);
        assert_eq!(result.y, y_res);
        assert_eq!(result.z, z_res);
    }

    #[test]
    fn test_complex_scalar_multiplication_and_sum_to_neutral() {
        let p: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
            16,
        )
        .unwrap();
        let a: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000fffffffffffffffffffffffc",
            16,
        )
        .unwrap();
        let b: BigInt = BigInt::from_str_radix(
            "5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
            16,
        )
        .unwrap();

        let curve: Curve = Curve::new(p, a, b).unwrap();

        let x1: BigInt = BigInt::from_str_radix(
            "33404865307973760982101572896420120260786884941743684464500659354489713797603",
            10,
        )
        .unwrap();
        let y1: BigInt = BigInt::from_str_radix(
            "36876760391286198122382864263676877270014017711839896656233225271501426068233",
            10,
        )
        .unwrap();
        let z1: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let scalar1: BigInt = BigInt::from_str_radix(
            "115792089210356248762697446949407573529996955224135760342422259061068512044367",
            10,
        )
        .unwrap();
        let scalar2: BigInt = BigInt::from_str_radix("2", 10).unwrap();

        let x_res: BigInt = BigInt::from_str_radix("0", 10).unwrap();
        let y_res: BigInt = BigInt::from_str_radix("1", 10).unwrap();
        let z_res: BigInt = BigInt::from_str_radix("0", 10).unwrap();

        let point1: ProjectivePoint = ProjectivePoint::new(x1.clone(), y1.clone(), z1.clone());

        let result: ProjectivePoint = curve.add(
            &curve.scalar_mul(&scalar1, &point1),
            &curve.scalar_mul(&scalar2, &point1),
        );

        assert_eq!(result.x, x_res);
        assert_eq!(result.y, y_res);
        assert_eq!(result.z, z_res);
    }

    #[test]
    fn test_complex_scalar_multiplication_and_sum_montgomery() {
        let p: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
            16,
        )
        .unwrap();
        let a: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000fffffffffffffffffffffffc",
            16,
        )
        .unwrap();
        let b: BigInt = BigInt::from_str_radix(
            "5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
            16,
        )
        .unwrap();

        let curve: Curve = Curve::new(p, a, b).unwrap();

        let x1: BigInt = BigInt::from_str_radix(
            "29424756536908666275014196764752012937549470486181629932252170285403408049517",
            10,
        )
        .unwrap();
        let y1: BigInt = BigInt::from_str_radix(
            "99699357898675284307950095025727949597877550412833207843633386200194164890561",
            10,
        )
        .unwrap();
        let z1: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let x2: BigInt = BigInt::from_str_radix(
            "95387701464026858343496773958108763791065521288270402154042613065295218611790",
            10,
        )
        .unwrap();
        let y2: BigInt = BigInt::from_str_radix(
            "56197753972877288816938393898610332375031967031226761314470230286071453671148",
            10,
        )
        .unwrap();
        let z2: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let scalar1: BigInt = BigInt::from_str_radix(
            "115166654073940935621181939619999211512712171442155919412146477633437056010824",
            10,
        )
        .unwrap();
        let scalar2: BigInt = BigInt::from_str_radix(
            "18076504642849865357726580190900882117966238578277710228703901228455155726225",
            10,
        )
        .unwrap();

        let x_res: BigInt = BigInt::from_str_radix(
            "37352616003680990833471240280964019197164191581612728497391994362995806189072",
            10,
        )
        .unwrap();
        let y_res: BigInt = BigInt::from_str_radix(
            "51702340847344877490410640050738219062093639452275966547398163436084404331575",
            10,
        )
        .unwrap();
        let z_res: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let point1: ProjectivePoint = ProjectivePoint::new(x1.clone(), y1.clone(), z1.clone());
        let point2: ProjectivePoint = ProjectivePoint::new(x2.clone(), y2.clone(), z2.clone());

        let computed_result_projective_point: ProjectivePoint = curve.add(
            &curve.scalar_mul_montgomery(&scalar1, &point1),
            &curve.scalar_mul_montgomery(&scalar2, &point2),
        );

        let computed_result_affine_point: AffinePoint =
            computed_result_projective_point.to_affine(&curve).unwrap();
        let result: ProjectivePoint = computed_result_affine_point.to_projective();

        assert_eq!(result.x, x_res);
        assert_eq!(result.y, y_res);
        assert_eq!(result.z, z_res);
    }

    #[test]
    fn test_complex_scalar_multiplication_and_sum_to_neutral_montgomery() {
        let p: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
            16,
        )
        .unwrap();
        let a: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000fffffffffffffffffffffffc",
            16,
        )
        .unwrap();
        let b: BigInt = BigInt::from_str_radix(
            "5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
            16,
        )
        .unwrap();

        let curve: Curve = Curve::new(p, a, b).unwrap();

        let x1: BigInt = BigInt::from_str_radix(
            "33404865307973760982101572896420120260786884941743684464500659354489713797603",
            10,
        )
        .unwrap();
        let y1: BigInt = BigInt::from_str_radix(
            "36876760391286198122382864263676877270014017711839896656233225271501426068233",
            10,
        )
        .unwrap();
        let z1: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let scalar1: BigInt = BigInt::from_str_radix(
            "115792089210356248762697446949407573529996955224135760342422259061068512044367",
            10,
        )
        .unwrap();
        let scalar2: BigInt = BigInt::from_str_radix("2", 10).unwrap();

        let x_res: BigInt = BigInt::from_str_radix("0", 10).unwrap();
        let y_res: BigInt = BigInt::from_str_radix("1", 10).unwrap();
        let z_res: BigInt = BigInt::from_str_radix("0", 10).unwrap();

        let point1: ProjectivePoint = ProjectivePoint::new(x1.clone(), y1.clone(), z1.clone());

        let result: ProjectivePoint = curve.add(
            &curve.scalar_mul_montgomery(&scalar1, &point1),
            &curve.scalar_mul_montgomery(&scalar2, &point1),
        );

        assert_eq!(result.x, x_res);
        assert_eq!(result.y, y_res);
        assert_eq!(result.z, z_res);
    }

    #[test]
    fn test_random_curve_mul_order() {
        let p: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
            16,
        )
        .unwrap();
        let a: BigInt = BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000fffffffffffffffffffffffc",
            16,
        )
        .unwrap();
        let b: BigInt = BigInt::from_str_radix(
            "5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
            16,
        )
        .unwrap();

        let curve: Curve = Curve::new(p, a, b).unwrap();

        let random_point: ProjectivePoint = curve.generate_random_projective_point_p256_naive();

        // println!("{:?}", random_point);

        let _ = &curve.scalar_mul_montgomery(&curve.get_order(), &random_point);
    }

    #[test]
    fn test_diffie_hellman_key_exchange() {
        let curve: Curve = helpers::get_curve_p256();
        let base_point: ProjectivePoint = curve.generate_random_projective_point_p256_fast();

        let alice: Actor = Actor::new(helpers::get_curve_p256());
        let bob: Actor = Actor::new(helpers::get_curve_p256());

        let alice_pre_key: ProjectivePoint = alice.generate_pre_key(&base_point);
        let bob_pre_key: ProjectivePoint = bob.generate_pre_key(&base_point);

        let alice_key: BigInt = alice.compute_common_secret(&bob_pre_key);
        let bob_key: BigInt = bob.compute_common_secret(&alice_pre_key);

        assert!(alice_key.eq(&bob_key));
    }
}

#[cfg(test)]
mod speed {
    use crate::{helpers, Actor, Curve, ProjectivePoint};
    use num_bigint::BigInt;
    use std::time::Instant;

    #[test]
    fn test_generate_random_projective_point() {
        let curve: Curve = helpers::get_curve_p256();

        let method_to_run = || -> u128 {
            let start_time: Instant = Instant::now();

            curve.generate_random_projective_point_p256_naive();

            let end_time: Instant = Instant::now();

            end_time.duration_since(start_time).as_nanos()
        };

        helpers::measure_average_execution_time(
            "generate_random_projective_point_p256",
            method_to_run,
            100,
        );
    }

    #[test]
    fn test_generate_random_projective_point_p256_fast() {
        let curve: Curve = helpers::get_curve_p256();

        let method_to_run = || -> u128 {
            let start_time: Instant = Instant::now();

            Curve::generate_random_projective_point_p256_fast(&curve);

            let end_time: Instant = Instant::now();

            end_time.duration_since(start_time).as_nanos()
        };

        helpers::measure_average_execution_time(
            "generate_random_projective_point_p256_fast",
            method_to_run,
            100,
        );
    }

    #[test]
    fn test_verify_affine_point() {
        let curve: Curve = helpers::get_curve_p256();

        let method_to_run = || -> u128 {
            let point: ProjectivePoint = curve.generate_random_projective_point_p256_naive();

            let start_time: Instant = Instant::now();

            curve.verify_affine_point(&point.to_affine(&curve).unwrap());

            let end_time: Instant = Instant::now();

            end_time.duration_since(start_time).as_nanos()
        };

        helpers::measure_average_execution_time("verify_affine_point", method_to_run, 100);
    }

    #[test]
    fn test_verify_projective_point() {
        let curve: Curve = helpers::get_curve_p256();

        let method_to_run = || -> u128 {
            let point: ProjectivePoint = curve.generate_random_projective_point_p256_naive();

            let start_time: Instant = Instant::now();

            curve.verify_projective_point(&point);

            let end_time: Instant = Instant::now();

            end_time.duration_since(start_time).as_nanos()
        };

        helpers::measure_average_execution_time("verify_projective_point", method_to_run, 100);
    }

    #[test]
    fn test_add() {
        let curve: Curve = helpers::get_curve_p256();

        let method_to_run = || -> u128 {
            let point1: ProjectivePoint = curve.generate_random_projective_point_p256_naive();
            let point2: ProjectivePoint = curve.generate_random_projective_point_p256_naive();

            let start_time: Instant = Instant::now();

            curve.add(&point1, &point2);

            let end_time: Instant = Instant::now();

            end_time.duration_since(start_time).as_nanos()
        };

        helpers::measure_average_execution_time("add", method_to_run, 100);
    }

    #[test]
    fn test_double() {
        let curve: Curve = helpers::get_curve_p256();

        let method_to_run = || -> u128 {
            let point: ProjectivePoint = curve.generate_random_projective_point_p256_naive();

            let start_time: Instant = Instant::now();

            curve.double(&point);

            let end_time: Instant = Instant::now();

            end_time.duration_since(start_time).as_nanos()
        };

        helpers::measure_average_execution_time("double", method_to_run, 100);
    }

    #[test]
    fn test_scalar_mul() {
        let curve: Curve = helpers::get_curve_p256();

        let method_to_run = || -> u128 {
            let point: ProjectivePoint = curve.generate_random_projective_point_p256_naive();
            let scalar: BigInt = helpers::generate_random_bigint(curve.p.bits() as usize);

            let start_time: Instant = Instant::now();

            curve.scalar_mul(&scalar, &point);

            let end_time: Instant = Instant::now();

            end_time.duration_since(start_time).as_nanos()
        };

        helpers::measure_average_execution_time("scalar_mul", method_to_run, 100);
    }

    #[test]
    fn test_scalar_mul_montgomery() {
        let curve: Curve = helpers::get_curve_p256();

        let method_to_run = || -> u128 {
            let point: ProjectivePoint = curve.generate_random_projective_point_p256_naive();
            let scalar: BigInt = helpers::generate_random_bigint(curve.p.bits() as usize);

            let start_time: Instant = Instant::now();

            curve.scalar_mul_montgomery(&scalar, &point);

            let end_time: Instant = Instant::now();

            end_time.duration_since(start_time).as_nanos()
        };

        helpers::measure_average_execution_time("scalar_mul_montgomery", method_to_run, 100);
    }

    #[test]
    fn test_diffie_hellman_pre_key_computation() {
        let method_to_run = || -> u128 {
            let alice: Actor = Actor::new(helpers::get_curve_p256());
            let curve: Curve = helpers::get_curve_p256();
            let base_point: ProjectivePoint = curve.generate_random_projective_point_p256_fast();

            let start_time: Instant = Instant::now();
            alice.generate_pre_key(&base_point);
            let end_time: Instant = Instant::now();

            end_time.duration_since(start_time).as_nanos()
        };

        helpers::measure_average_execution_time("diffie_hellman_pre_key_computation", method_to_run, 100);
    }

    #[test]
    fn test_diffie_hellman_key_computation() {
        let method_to_run = || -> u128 {
            let alice: Actor = Actor::new(helpers::get_curve_p256());
            let curve: Curve = helpers::get_curve_p256();
            let bob_pre_key: ProjectivePoint = curve.generate_random_projective_point_p256_fast();

            let start_time: Instant = Instant::now();
            alice.compute_common_secret(&bob_pre_key);
            let end_time: Instant = Instant::now();

            end_time.duration_since(start_time).as_nanos()
        };

        helpers::measure_average_execution_time("diffie_hellman_key_computation", method_to_run, 100);
    }
}
