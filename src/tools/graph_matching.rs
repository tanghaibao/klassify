use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};

/// Greedy disjoint maximum-*weight* selection (fast, not exact).
/// Keeps pairs disjoint by excluding nodes once used.
pub fn greedy_max_weight_matching(
    pair_to_reads: &HashMap<(String, String), Vec<String>>,
    min_read_support: usize,
) -> Vec<(String, String)> {
    #[derive(Clone)]
    struct Edge {
        a: String,
        b: String,
        w: usize,
    }
    let mut edges: Vec<Edge> = pair_to_reads
        .iter()
        .map(|((a, b), reads)| Edge {
            a: a.clone(),
            b: b.clone(),
            w: reads.len(),
        })
        .collect();

    // sort by weight desc, then tie-break lexicographically for stability
    edges.sort_by(|x, y| match y.w.cmp(&x.w) {
        Ordering::Equal => (x.a.as_str(), x.b.as_str()).cmp(&(y.a.as_str(), y.b.as_str())),
        other => other,
    });

    let mut used: HashSet<String> = HashSet::new();
    let mut out: Vec<(String, String)> = Vec::new();

    for e in edges {
        if e.w < min_read_support {
            continue;
        }
        if used.contains(&e.a) || used.contains(&e.b) {
            continue;
        }
        used.insert(e.a.clone());
        used.insert(e.b.clone());
        // canonicalize ordering
        if e.a <= e.b {
            out.push((e.a, e.b));
        } else {
            out.push((e.b, e.a));
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_greedy_basic() {
        let mut p2r: HashMap<(String, String), Vec<String>> = HashMap::new();
        // A-B:3, A-C:1, D-C:2  -> greedy chooses (A,B) and (C,D)
        p2r.insert(
            ("A:1-2".into(), "B:3-4".into()),
            vec!["r1".into(), "r2".into(), "r3".into()],
        );
        p2r.insert(("A:1-2".into(), "C:5-6".into()), vec!["r4".into()]);
        p2r.insert(
            ("D:7-8".into(), "C:5-6".into()),
            vec!["r5".into(), "r6".into()],
        );

        let sel = greedy_max_weight_matching(&p2r, 1);
        let set: HashSet<_> = sel.into_iter().collect();
        assert!(set.contains(&("A:1-2".to_string(), "B:3-4".to_string())));
        assert!(set.contains(&("C:5-6".to_string(), "D:7-8".to_string())));
        assert_eq!(set.len(), 2);
    }
}
