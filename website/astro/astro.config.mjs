import { defineConfig } from 'astro/config';
import tailwindcss from '@tailwindcss/vite';
import sitemap from '@astrojs/sitemap';

// The site is served from the GitHub Pages custom domain at the root.
// (galaxyproject.github.io/tools-iuc redirects to it, stripping the prefix.)
export default defineConfig({
  site: 'https://iuc.galaxyproject.org',
  base: '/',
  output: 'static',
  integrations: [sitemap()],
  trailingSlash: 'ignore',
  vite: {
    plugins: [tailwindcss()],
  },
});
