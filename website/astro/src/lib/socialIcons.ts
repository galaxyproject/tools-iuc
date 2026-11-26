// Link builders for contributor social profiles.
//
// Brand icon artwork lives in src/components/mdx/brand-icons.ts — a verbatim
// copy of galaxy-hub's astro/src/components/mdx/brand-icons.ts. To update the
// icons, re-copy that file instead of editing it here. `localIcons` holds the
// few non-brand marks used alongside the brand icons (same full-SVG string
// format, so BrandIcon.astro can look both up by name).

// Filled "website" mark — not a brand icon; galaxy-hub uses lucide's globe
// for these links instead.
export const localIcons: Record<string, string> = {
  Globe: `<svg
  xmlns="http://www.w3.org/2000/svg"
  width="24"
  height="24"
  viewBox="0 0 24 24"
  fill="currentColor"
>
  <path d="M12 0C5.373 0 0 5.373 0 12s5.373 12 12 12 12-5.373 12-12S18.627 0 12 0zm0 22c-5.514 0-10-4.486-10-10S6.486 2 12 2s10 4.486 10 10-4.486 10-10 10zm-1-19.5C8.1 3.4 5.5 7 5.5 12s2.6 8.6 5.5 9.5V2.5zm2 0v19c2.9-.9 5.5-4.5 5.5-9.5S15.9 3.4 13 2.5z" />
</svg>`,
};

export function socialLink(type: string, value: string): string | null {
  switch (type) {
    case "github":
      return `https://github.com/${value}`;
    case "orcid":
      return value.startsWith("http") ? value : `https://orcid.org/${value}`;
    case "linkedin":
      return `https://www.linkedin.com/in/${value}`;
    case "matrix":
      return `https://matrix.to/#/${value}`;
    case "fediverse":
    case "mastodon":
      return value.startsWith("http") ? value : `https://${value}`;
    case "bluesky":
      return `https://bsky.app/profile/${value}`;
    case "twitter":
      return `https://twitter.com/${value}`;
    case "googleScholar":
    case "google-scholar":
      return `https://scholar.google.com/citations?user=${value}`;
    case "researchgate":
      return `https://www.researchgate.net/profile/${value}`;
    case "globe":
    case "website":
    case "url":
      return value;
    default:
      return null;
  }
}
